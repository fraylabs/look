from OCP.BRepPrimAPI import BRepPrimAPI_MakeCylinder,BRepPrimAPI_MakeCone,BRepPrimAPI_MakeBox
from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse,BRepAlgoAPI_Cut
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCP.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCP.gp import gp_Ax1,gp_Ax2,gp_Pnt,gp_Dir,gp_Trsf
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_EDGE
from OCP.TopoDS import TopoDS
from OCP.STEPControl import STEPControl_Writer,STEPControl_AsIs
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.IFSelect import IFSelect_RetDone
from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.GeomAbs import GeomAbs_Circle
import pathlib,math,sys
out=pathlib.Path(__file__).resolve().parent;out.mkdir(exist_ok=True)
cylinder=BRepPrimAPI_MakeCylinder(5,8).Shape()
base=BRepPrimAPI_MakeCylinder(5,4).Shape();cap=BRepPrimAPI_MakeCone(gp_Ax2(gp_Pnt(0,0,4),gp_Dir(0,0,1)),5,3,3).Shape();fuse=BRepAlgoAPI_Fuse(base,cap);fuse.Build();shared=fuse.Shape()
box=BRepPrimAPI_MakeBox(12,8,6).Shape();fillet=BRepFilletAPI_MakeFillet(box);explorer=TopExp_Explorer(box,TopAbs_EDGE)
while explorer.More():fillet.Add(0.8,TopoDS.Edge_s(explorer.Current()));explorer.Next()
fillet.Build();assert fillet.IsDone()
outer=BRepPrimAPI_MakeCylinder(1,1).Shape();inner=BRepPrimAPI_MakeCylinder(0.9995,1).Shape()
rotation=gp_Trsf();rotation.SetRotation(gp_Ax1(gp_Pnt(0,0,0),gp_Dir(0,0,1)),math.pi/24)
inner=BRepBuilderAPI_Transform(inner,rotation,True).Shape()
annulus=BRepAlgoAPI_Cut(outer,inner);annulus.Build();assert annulus.IsDone()
long_cylinder=BRepPrimAPI_MakeCylinder(4.25,1000).Shape();rounded=BRepFilletAPI_MakeFillet(long_cylinder);explorer=TopExp_Explorer(long_cylinder,TopAbs_EDGE)
while explorer.More():
 edge=TopoDS.Edge_s(explorer.Current())
 if BRepAdaptor_Curve(edge).GetType()==GeomAbs_Circle:rounded.Add(0.5,edge)
 explorer.Next()
rounded.Build();assert rounded.IsDone()
rounded_shape=BRepBuilderAPI_Transform(rounded.Shape(),rotation,True).Shape()
for name,shape in [('cylinder-seam',cylinder),('shared-curved-edge',shared),('filleted-box',fillet.Shape()),('thin-annulus',annulus.Shape()),('torus-fillet',rounded_shape)]:
 if len(sys.argv)>1 and name not in sys.argv[1:]:continue
 assert BRepCheck_Analyzer(shape).IsValid(),name
 writer=STEPControl_Writer();writer.Transfer(shape,STEPControl_AsIs);assert writer.Write(str(out/(name+'.step')))==IFSelect_RetDone
 print(name,flush=True)
