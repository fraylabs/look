These four synthetic solids were authored for the mesh conformity regressions
using OCCT 7.9 via OCP. They contain no third-party CAD geometry and are licensed
under the repository's MIT OR Apache-2.0 terms. `generate.py` records the
construction; all four passed `BRepCheck_Analyzer.IsValid()` before STEP export.

- `cylinder-seam.step`: radius 5, height 8, with a periodic surface seam.
- `shared-curved-edge.step`: radius 5, height 4 cylinder fused to a cone of height
  3 whose radii are 5 and 3. Their curved circular edge is shared.
- `filleted-box.step`: 12 × 8 × 6 box, with radius 0.8 on all twelve edges.
- `thin-annulus.step`: radius 1 outer cylinder minus a radius 0.9995 inner
  cylinder, height 1. The inner seam is rotated by 7.5 degrees, so coarse
  polygonal boundary approximations cross despite the valid curved boundaries.
