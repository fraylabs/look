//! External-corpus validation harness. No source CAD is embedded here.
use serde_json::{Value, json};
use truck_meshalgo::tessellation::MeshableShape;
use truck_stepio::r#in::{Table, step_geometry::*};

fn point(value: &Value) -> Point3 {
    Point3::new(
        value[0].as_f64().unwrap(),
        value[1].as_f64().unwrap(),
        value[2].as_f64().unwrap(),
    )
}
fn main() {
    let path = std::env::args().nth(1).unwrap();
    let input: Vec<Value> = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let mut output = Vec::new();
    for row in input {
        let table =
            Table::from_step(&std::fs::read_to_string(row["path"].as_str().unwrap()).unwrap())
                .unwrap();
        let sid = row["shell"].as_u64().unwrap();
        let (shell, losses) = table
            .to_compressed_shell_with_losses(sid, &table.shell[&sid])
            .unwrap();
        let surface = &shell.faces[0].surface;
        let Surface::ElementarySurface(ElementarySurface::ToroidalSurface(processor)) = &surface
        else {
            panic!("expected torus")
        };
        let swapped = !processor.orientation();
        let samples = row["samples"]
            .as_array()
            .unwrap()
            .iter()
            .map(|sample| {
                let u = sample["uv"][0].as_f64().unwrap();
                let v = sample["uv"][1].as_f64().unwrap();
                let uv = if swapped { (v, u) } else { (u, v) };
                let expected = point(&sample["point"]);
                let inverse = surface.search_parameter(expected, None, 100);
                let nearest = surface.search_nearest_parameter(expected, None, 100);
                json!({"point_error":surface.subs(uv.0,uv.1).distance(expected),
                "normal":surface.normal(uv.0,uv.1).normalize(),
                "inverse_error":inverse.map(|(u,v)|surface.subs(u,v).distance(expected)),
                "nearest_error":nearest.map(|(u,v)|surface.subs(u,v).distance(expected))})
            })
            .collect::<Vec<_>>();
        let boundary = row["boundary"].as_array().unwrap().iter().map(|p| {
            let expected=point(p);
            json!({"inverse_error":surface.search_parameter(expected,None,100).map(|(u,v)|surface.subs(u,v).distance(expected)),
                "nearest_error":surface.search_nearest_parameter(expected,None,100).map(|(u,v)|surface.subs(u,v).distance(expected))})
        }).collect::<Vec<_>>();
        let mesh = shell.triangulation(0.001);
        let polygon = mesh.faces.first().and_then(|f| f.surface.as_ref());
        let mesh_data = polygon.map(|poly| {
            let positions = poly.positions();
            let triangles = poly
                .tri_faces()
                .iter()
                .map(|f| f.map(|v| v.pos))
                .collect::<Vec<_>>();
            let area = triangles
                .iter()
                .map(|f| {
                    (positions[f[1]] - positions[f[0]])
                        .cross(positions[f[2]] - positions[f[0]])
                        .magnitude()
                        / 2.0
                })
                .sum::<f64>();
            json!({"positions":positions,"triangles":triangles,"area":area})
        });
        output.push(
            json!({"model":row["model"],"face":row["face"],"surface":row["surface"],
            "samples":samples,"boundary":boundary,"losses":format!("{losses:?}"),"mesh":mesh_data}),
        );
    }
    println!("{}", serde_json::to_string(&output).unwrap());
}
