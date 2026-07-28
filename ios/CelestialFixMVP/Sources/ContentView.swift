// ContentView.swift — minimal MVP UI: tele preview, Capture, and the fix + covariance.
// © 2026 MIT.

import SwiftUI
import AVFoundation
import MapKit

/// UIKit preview layer wrapped for SwiftUI.
struct CameraPreview: UIViewRepresentable {
    let session: AVCaptureSession
    func makeUIView(context: Context) -> PreviewView {
        let v = PreviewView(); v.videoPreviewLayer.session = session
        v.videoPreviewLayer.videoGravity = .resizeAspectFill
        return v
    }
    func updateUIView(_ uiView: PreviewView, context: Context) {}
    final class PreviewView: UIView {
        override class var layerClass: AnyClass { AVCaptureVideoPreviewLayer.self }
        var videoPreviewLayer: AVCaptureVideoPreviewLayer { layer as! AVCaptureVideoPreviewLayer }
    }
}

struct ContentView: View {
    @StateObject private var capture: CaptureController
    @StateObject private var vm: FixViewModel
    @State private var region = MKCoordinateRegion(
        center: .init(latitude: 51.5, longitude: 0),
        span: .init(latitudeDelta: 2, longitudeDelta: 2))

    init() {
        let c = CaptureController()
        _capture = StateObject(wrappedValue: c)
        _vm = StateObject(wrappedValue: FixViewModel(capture: c))
    }

    var body: some View {
        VStack(spacing: 0) {
            ZStack(alignment: .center) {
                CameraPreview(session: capture.session).ignoresSafeArea(edges: .top)
                // crosshair on the principal point
                Crosshair().stroke(.white.opacity(0.6), lineWidth: 1).frame(width: 40, height: 40)
            }
            .frame(maxHeight: .infinity)

            resultPanel
                .padding()
                .background(.ultraThinMaterial)
        }
        .onAppear { capture.start() }
        .onDisappear { capture.stop() }
    }

    @ViewBuilder private var resultPanel: some View {
        VStack(alignment: .leading, spacing: 10) {
            Text(vm.status).font(.callout).foregroundStyle(.primary)
            if let r = vm.result {
                Map(coordinateRegion: $region,
                    annotationItems: [FixPin(coord: .init(latitude: r.fix.latDeg,
                                                          longitude: r.fix.lonDeg))]) { pin in
                    MapMarker(coordinate: pin.coord, tint: .cyan)
                }
                .frame(height: 140).clipShape(RoundedRectangle(cornerRadius: 10))
                Grid(alignment: .leading, horizontalSpacing: 14, verticalSpacing: 3) {
                    row("Moon altitude", String(format: "%.3f° ± %.3f°", r.measAltDeg, r.sigmaAltDeg))
                    row("Plate scale", String(format: "%.2f″/px  (R=%.0f px)", r.arcsecPerPixel, r.radiusPx))
                    row("Limb fit RMSE", String(format: "%.3f px", r.rmsePx))
                    row("Fix 1σ (major)", String(format: "%.1f km", r.fix.oneSigmaMajorKm))
                    row("GN iterations", "\(r.fix.iterations)")
                }.font(.system(.footnote, design: .monospaced))
            }
            Button {
                vm.onCapture()
                if let r = vm.result {
                    region.center = .init(latitude: r.fix.latDeg, longitude: r.fix.lonDeg)
                }
            } label: {
                Label("Capture sight", systemImage: "camera.aperture")
                    .font(.headline).frame(maxWidth: .infinity)
            }
            .buttonStyle(.borderedProminent)
        }
    }

    private func row(_ k: String, _ v: String) -> some View {
        GridRow { Text(k).foregroundStyle(.secondary); Text(v) }
    }
}

struct FixPin: Identifiable { let id = UUID(); let coord: CLLocationCoordinate2D }

struct Crosshair: Shape {
    func path(in r: CGRect) -> Path {
        var p = Path()
        p.move(to: CGPoint(x: r.midX, y: r.minY)); p.addLine(to: CGPoint(x: r.midX, y: r.maxY))
        p.move(to: CGPoint(x: r.minX, y: r.midY)); p.addLine(to: CGPoint(x: r.maxX, y: r.midY))
        return p
    }
}
