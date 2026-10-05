use std::sync::atomic::{AtomicUsize, Ordering};
use truck_geometry::prelude::*;
use truck_meshalgo::prelude::*;
use truck_meshalgo::tessellation::{formal, TessellationFailureReason};
use truck_modeling::{Curve, Surface};
use truck_topology::compress::*;

fn rectangle(height: f64) -> CompressedShell<Point3, Curve, Surface> {
    let vertices = vec![
        Point3::origin(),
        Point3::new(1., 0., 0.),
        Point3::new(1., height, 0.),
        Point3::new(0., height, 0.),
    ];
    let edges = (0..4)
        .map(|i| CompressedEdge {
            vertices: (i, (i + 1) % 4),
            curve: Curve::Line(Line(vertices[i], vertices[(i + 1) % 4])),
        })
        .collect();
    CompressedShell {
        vertices,
        edges,
        source_geometric_uncertainty: None,
        faces: vec![CompressedFace {
            boundaries: vec![(0..4)
                .map(|index| CompressedEdgeIndex {
                    index,
                    orientation: true,
                })
                .collect()],
            orientation: true,
            provenance: Default::default(),
            surface: Surface::Plane(Plane::new(
                Point3::origin(),
                Point3::new(1., 0., 0.),
                Point3::new(0., 1., 0.),
            )),
        }],
    }
}

fn with_retry(
    shell: &CompressedShell<Point3, Curve, Surface>,
    calls: &AtomicUsize,
) -> truck_meshalgo::tessellation::MeshedShellOutcome {
    shell.robust_triangulation_with_inverse_retry_outcome(
        0.01,
        unevidenced_lattice,
        |_| {
            formal::SupportSurfaceSchema::not_structurally_identified(
                formal::SchemaIdentificationFailure::NoStructuralReader {
                    representation: "inverse_retry_test",
                },
            )
        },
        |_| {
            formal::CurveSchema::not_structurally_identified(
                formal::CurveSchemaFailure::NoStructuralReader {
                    representation: "inverse_retry_test",
                },
            )
        },
        |_| Err("no_cylinder"),
        |_| {
            formal::CurveSchema::not_structurally_identified(
                formal::CurveSchemaFailure::NoStructuralReader {
                    representation: "inverse_retry_test",
                },
            )
        },
        |_| None,
        |_| Err("no_cone"),
        |_| Err("no_torus"),
        |surface| {
            calls.fetch_add(1, Ordering::SeqCst);
            Some(surface.clone())
        },
    )
}

#[test]
fn inverse_retry_never_touches_a_successful_face_or_its_shared_edges() {
    let shell = rectangle(1.);
    let baseline = shell.robust_triangulation_with_lattice_outcome(0.01, unevidenced_lattice);
    let calls = AtomicUsize::new(0);
    let actual = with_retry(&shell, &calls);
    assert_eq!(calls.load(Ordering::SeqCst), 0);
    assert!(actual.face_failures[0].is_none());
    assert!(actual.shell.faces[0]
        .surface
        .as_ref()
        .is_some_and(|mesh| !mesh.tri_faces().is_empty()));
    assert_eq!(
        actual.shell.faces[0].surface,
        baseline.shell.faces[0].surface
    );
    for (a, b) in actual.shell.edges.iter().zip(&baseline.shell.edges) {
        assert_eq!(a.vertices, b.vertices);
        assert_eq!(a.curve, b.curve);
    }
}

#[test]
fn inverse_retry_cannot_recover_a_certified_zero_area_trim() {
    let shell = rectangle(0.);
    let calls = AtomicUsize::new(0);
    let actual = with_retry(&shell, &calls);
    assert_eq!(calls.load(Ordering::SeqCst), 0);
    assert_eq!(
        actual.face_failures[0].as_ref().map(|f| f.reason),
        Some(TessellationFailureReason::RejectedDegenerate)
    );
}
