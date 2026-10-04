from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder,BRepPrimAPI_MakeCone,BRepPrimAPI_MakeBox
from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCP.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCP.gp import gp_Ax2,gp_Pnt,gp_Dir
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_EDGE
from OCP.TopoDS import TopoDS
from OCP.STEPControl import STEPControl_Writer,STEPControl_AsIs
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.IFSelect import IFSelect_RetDone
import pathlib
out=pathlib.Path(__file__).resolve().parent;out.mkdir(exist_ok=True)
cylinder=BRepPrimAPI_MakeCylinder(5,8).Shape()
base=BRepPrimAPI_MakeCylinder(5,4).Shape();cap=BRepPrimAPI_MakeCone(gp_Ax2(gp_Pnt(0,0,4),gp_Dir(0,0,1)),5,3,3).Shape();fuse=BRepAlgoAPI_Fuse(base,cap);fuse.Build();shared=fuse.Shape()
box=BRepPrimAPI_MakeBox(12,8,6).Shape();fillet=BRepFilletAPI_MakeFillet(box);explorer=TopExp_Explorer(box,TopAbs_EDGE)
while explorer.More():fillet.Add(0.8,TopoDS.Edge_s(explorer.Current()));explorer.Next()
fillet.Build();assert fillet.IsDone()
for name,shape in [('cylinder-seam',cylinder),('shared-curved-edge',shared),('filleted-box',fillet.Shape())]:
 assert BRepCheck_Analyzer(shape).IsValid(),name
 writer=STEPControl_Writer();writer.Transfer(shape,STEPControl_AsIs);assert writer.Write(str(out/(name+'.step')))==IFSelect_RetDone
 print(name,flush=True)
