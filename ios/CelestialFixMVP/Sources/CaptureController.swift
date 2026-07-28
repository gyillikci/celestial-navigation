// CaptureController.swift — tele-camera capture + IMU gravity at the shutter instant.
// Grabs one telephoto frame as luminance and the CoreMotion gravity vector in the
// device frame, so the pipeline can turn them into an altitude sight.
//
// MVP: single-camera (telephoto), single frame.  The full app adds the ultrawide
// (optical horizon), the least-rotation gate, and multi-cam.  © 2026 MIT.

import AVFoundation
import CoreMotion
import CoreVideo
import simd

public final class CaptureController: NSObject, ObservableObject,
                                      AVCaptureVideoDataOutputSampleBufferDelegate {

    public let session = AVCaptureSession()
    private let output = AVCaptureVideoDataOutput()
    private let queue = DispatchQueue(label: "capture.video")
    private let motion = CMMotionManager()

    private var latest: CVPixelBuffer?
    private let bufLock = NSLock()

    @Published public var ready = false
    @Published public var errorText: String?

    public override init() {
        super.init()
        configure()
        startMotion()
    }

    private func configure() {
        session.beginConfiguration()
        session.sessionPreset = .photo
        guard let device = AVCaptureDevice.default(.builtInTelephotoCamera,
                                                   for: .video, position: .back)
                ?? AVCaptureDevice.default(.builtInWideAngleCamera, for: .video, position: .back),
              let input = try? AVCaptureDeviceInput(device: device),
              session.canAddInput(input) else {
            errorText = "No telephoto camera available."
            session.commitConfiguration(); return
        }
        session.addInput(input)
        output.videoSettings = [kCVPixelBufferPixelFormatTypeKey as String:
                                kCVPixelFormatType_32BGRA]
        output.setSampleBufferDelegate(self, queue: queue)
        if session.canAddOutput(output) { session.addOutput(output) }
        // Lock exposure/focus so the disk is stable (user sets it on the Moon).
        try? device.lockForConfiguration()
        if device.isFocusModeSupported(.locked) { device.focusMode = .locked }
        device.unlockForConfiguration()
        session.commitConfiguration()
        ready = true
    }

    public func start() { if !session.isRunning { queue.async { self.session.startRunning() } } }
    public func stop()  { if session.isRunning { session.stopRunning() } }

    private func startMotion() {
        guard motion.isDeviceMotionAvailable else { return }
        motion.deviceMotionUpdateInterval = 1.0 / 100.0
        motion.startDeviceMotionUpdates(using: .xArbitraryZVertical)
    }

    /// Gravity unit vector in the DEVICE frame at this instant (points to ground).
    public var gravity: SIMD3<Double>? {
        guard let g = motion.deviceMotion?.gravity else { return nil }
        return simd_normalize(SIMD3<Double>(g.x, g.y, g.z))
    }

    public func captureLuminance(completion: @escaping (GrayImage?, SIMD3<Double>?) -> Void) {
        bufLock.lock(); let pb = latest; bufLock.unlock()
        let g = gravity
        guard let pb else { completion(nil, g); return }
        completion(Self.toGray(pb), g)
    }

    public func captureOutput(_ o: AVCaptureOutput, didOutput sampleBuffer: CMSampleBuffer,
                              from connection: AVCaptureConnection) {
        guard let pb = CMSampleBufferGetImageBuffer(sampleBuffer) else { return }
        bufLock.lock(); latest = pb; bufLock.unlock()
    }

    /// BGRA CVPixelBuffer → luminance GrayImage (Rec.601).  Optionally strides for speed.
    static func toGray(_ pb: CVPixelBuffer, stride s: Int = 1) -> GrayImage? {
        CVPixelBufferLockBaseAddress(pb, .readOnly)
        defer { CVPixelBufferUnlockBaseAddress(pb, .readOnly) }
        guard let base = CVPixelBufferGetBaseAddress(pb) else { return nil }
        let w = CVPixelBufferGetWidth(pb), h = CVPixelBufferGetHeight(pb)
        let rowBytes = CVPixelBufferGetBytesPerRow(pb)
        let ptr = base.assumingMemoryBound(to: UInt8.self)
        let ow = w / s, oh = h / s
        var out = [Float](repeating: 0, count: ow * oh)
        for oy in 0..<oh {
            let row = ptr + (oy * s) * rowBytes
            for ox in 0..<ow {
                let p = row + (ox * s) * 4          // BGRA
                let b = Float(p[0]), g = Float(p[1]), r = Float(p[2])
                out[oy * ow + ox] = 0.114*b + 0.587*g + 0.299*r
            }
        }
        return GrayImage(width: ow, height: oh, pixels: out)
    }
}
