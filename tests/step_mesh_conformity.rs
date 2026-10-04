//! Synthetic OCCT-valid solids must retain closed, oriented mesh boundaries.
use look::{config::UpAxis, scene::compile_scene, timing::Timings};
use std::collections::{HashMap, HashSet};

fn assert_closed_fixture(name: &str) {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/mesh-conformity")
        .join(format!("{name}.step"));
    let scene = compile_scene(&path, UpAxis::Z, &mut Timings::default()).unwrap();
    assert!(scene.assembly_structure_errors.is_empty());
    assert_eq!(scene.statistics.step_import.unwrap().lost_faces, 0);
    assert_eq!(scene.geometries.len(), 1);
    let geometry = &scene.geometries[0];
    // The fixtures span millimetres near the origin. This is a numerical weld
    // below the STEP source uncertainty, independent of the meshing tolerance.
    let key = |position: [f32; 3]| position.map(|value| (f64::from(value) * 1e6).round() as i64);
    let mut edges: HashMap<[[i64; 3]; 2], (usize, i32)> = HashMap::new();
    let mut triangles = HashSet::new();
    for triangle in geometry.indices.chunks_exact(3) {
        let points = [triangle[0], triangle[1], triangle[2]]
            .map(|index| key(geometry.vertices[index as usize].position));
        if points[0] == points[1] || points[1] == points[2] || points[2] == points[0] {
            continue;
        }
        let mut ordered = points;
        ordered.sort_unstable();
        assert!(triangles.insert(ordered), "{name}: duplicate triangle");
        for [a, b] in [
            [points[0], points[1]],
            [points[1], points[2]],
            [points[2], points[0]],
        ] {
            let (edge, direction) = if a < b { ([a, b], 1) } else { ([b, a], -1) };
            let count = edges.entry(edge).or_default();
            count.0 += 1;
            count.1 += direction;
        }
    }
    assert!(!triangles.is_empty());
    let open = edges.values().filter(|count| count.0 == 1).count();
    let nonmanifold = edges.values().filter(|count| count.0 > 2).count();
    let inconsistent = edges.values().filter(|count| count.1 != 0).count();
    assert_eq!((open, nonmanifold, inconsistent), (0, 0, 0), "{name}");
}

#[test]
fn cylinder_seam_is_closed_and_manifold() {
    assert_closed_fixture("cylinder-seam");
}

#[test]
fn two_faces_sharing_a_curved_edge_are_closed_and_manifold() {
    assert_closed_fixture("shared-curved-edge");
}

#[test]
fn filleted_box_is_closed_and_manifold() {
    assert_closed_fixture("filleted-box");
}

#[test]
fn phase_offset_thin_annulus_is_closed_and_manifold() {
    assert_closed_fixture("thin-annulus");
}

#[test]
fn circular_edges_shared_with_tori_keep_the_angular_floor() {
    assert_closed_fixture("torus-fillet");
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/mesh-conformity/torus-fillet.step");
    let scene = compile_scene(&path, UpAxis::Z, &mut Timings::default()).unwrap();
    let points = &scene.geometries[0].vertices;
    // The analytically authored carrier has radius 4.25. The source is rotated
    // by half the angular floor's interval; all circular joins meet tori.
    for axis in [0, 2] {
        let extent = points
            .iter()
            .map(|v| v.position[axis].abs())
            .fold(0.0_f32, f32::max);
        assert!(4.25 - extent < 0.05, "axis {axis}: extent {extent}");
    }
}
