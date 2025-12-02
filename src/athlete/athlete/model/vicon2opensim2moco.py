import numpy as np
import opensim

from opensim import Model

opensim.ModelVisualizer.addDirToGeometrySearchPaths("/home/lisca/biomec/sim/src/opensim-models/Geometry")

class Vicon2OpenSim2MocoModelConverter:

    _map_markers = {

        "L.Toe"   : "LTOE",
        "R.Toe"   : "RTOE",

        "L.Ankle" : "LANK",
        "R.Ankle" : "RANK",

        "L.Heel"  : "LHEE",
        "R.Heel"  : "RHEE",

        "L.SH1"   : "LTIB",
        "R.SH1"   : "RTIB",
        "L.SH2"   : None,
        "R.SH2"   : None,
        "L.SH3"   : None,
        "R.SH3"   : None,
        "R.SH4"   : None,

        "L.Knee"  : "LKNE",
        "R.Knee"  : "RKNE",

        "L.TH1"   : "LTHI",
        "R.TH1"   : "RTHI",
        "L.TH2"   : None,
        "R.TH2"   : None,
        "L.TH3"   : None,
        "R.TH3"   : None,
        "L.TH4"   : None,

        "L.ASIS"  : "LASI",
        "R.ASIS"  : "RASI",

        "L.PSIS"  : "LPSI",
        "R.PSIS"  : "RPSI",

        "L.MT5"   : None,
        "R.MT5"   : None,

        "S2"      : None,

    }

    _file_model_opensim_matlab_tools = None
    _file_model_moco                 = None
    _file_model_vicon2opensim2moco   = None

    _model_opensim_matlab_tools      = None
    _model_moco                      = None
    _model_vicon2opensim2moco        = None

    def __init__(self,
        file_model_opensim_matlab_tools,
        file_model_moco,
        file_model_vicon2opensim2moco):

        self._file_model_opensim_ematlab_tools = file_model_opensim_matlab_tools
        self._file_model_moco                  = file_model_moco
        # Save the modified copy of the Moco model into a different file.
        self._file_model_vicon2opensim2moco    = file_model_vicon2opensim2moco

        self._model_opensim_matlab_tools = opensim.Model(self._file_model_opensim_ematlab_tools)
        self._model_moco                 = opensim.Model(self._file_model_moco)
        # Modify a copy of the Moco model.
        self._model_vicon2opensim2moco   = opensim.Model(self._file_model_moco)

    def adjust_marker_used(self,):

        marker_set = self._model_vicon2opensim2moco.updMarkerSet()
        marker_removed = []
        for marker in marker_set:
            name_old = marker.getName()
            name_new = self._map_markers[name_old]
            if name_new is not None:
                print(f"{name_old:>8} -> {name_new}")
                marker.setName(name_new)
            else:
                print(f"{name_old:>8} -> removed")
                marker_removed.append(marker)

        print("Markers removed:")
        for marker in marker_removed:
            print(f"{marker.getName():>8}")
            marker_set.remove(marker)

        print("Markers remaining: ")
        for marker in marker_set:
            print(f"{marker.getName():>8}")

    def adjust_marker_positions(self):

        marker_set = self._model_vicon2opensim2moco.updMarkerSet()

        marker_set.getComponent("LTHI").set_location(opensim.Vec3( 0.06955, -0.23318, -0.06940))
        marker_set.getComponent("RTHI").set_location(opensim.Vec3( 0.06909, -0.23350,  0.06997))

        marker_set.getComponent("LTIB").set_location(opensim.Vec3( 0.03904, -0.25787, -0.05331))
        marker_set.getComponent("RTIB").set_location(opensim.Vec3( 0.03979, -0.25752,  0.05310))

        marker_set.getComponent("LANK").set_location(opensim.Vec3(-0.01256, -0.45530, -0.05330))
        marker_set.getComponent("RANK").set_location(opensim.Vec3(-0.01271, -0.45539,  0.05346))

        marker_set.getComponent("LTOE").set_location(opensim.Vec3( 0.22775,  0.02519,  0.02538))
        marker_set.getComponent("RTOE").set_location(opensim.Vec3( 0.22776,  0.02518, -0.02593))

    def adjust_joints(self):

#        joint_set_moco  = self._model_moco.updJointSet()
#        joint_set_omt   = self._model_opensim_matlab_tools.updJointSet()
#        joint_set_v2o2m = self._model_vicon2opensim2moco.updJointSet()
#
#        print("Swapping in the joints of Moco's model:")
#
#        joint_set_v2o2m.remove(13)
#        joint_set_v2o2m.insert(13, joint_set_omt.get(9))
#
#        joint_set_v2o2m.remove(12)
#        joint_set_v2o2m.insert(12, joint_set_omt.get(4))
#
#        joint_set_v2o2m.remove(9)
#        joint_set_v2o2m.insert(9, joint_set_omt.get(10))
#
#        joint_set_v2o2m.remove(10)
#        joint_set_v2o2m.insert(10, joint_set_omt.get(5))

#        # subtalar_r_v2o2m <- subtalar_r_omt
#        print(f"{joint_set_v2o2m.get(12).getName()} <- {joint_set_omt.get(4).getName()}")
#        joint = joint_set_v2o2m.get(12)
#        joint = joint_set_omt.get(4)
#
#        # subtalar_l_v2o2m <- subtalar_l_omt
#        print(f"{joint_set_v2o2m.get(13).getName()} <- {joint_set_omt.get(9).getName()}")
#        joint = joint_set_v2o2m.get(13)
#        joint = joint_set_omt.get(9)
#
#        # mtp_r_v2o2m <- mtp_r_omt
#        print(f"{joint_set_v2o2m.get(10).getName()} <- {joint_set_omt.get(5).getName()}")
#        joint = joint_set_v2o2m.get(10)
#        joint = joint_set_omt.get(5)
#
#        # mpt_l_v2o2m <- mpt_l_omt
#        print(f"{joint_set_v2o2m.get(9).getName()} <- {joint_set_omt.get(10).getName()}")
#        joint = joint_set_v2o2m.get(9)
#        joint = joint_set_omt.get(10)
#
#        self._model_vicon2opensim2moco.finalizeConnections();

        pass

    def save_to_file(self,):
        self._model_vicon2opensim2moco.printToXML(self._file_model_vicon2opensim2moco)

if __name__ == "__main__":

    file_model_opensim_matlab_tools = "/home/lisca/biomec/sim/src/athlete/matlab/vicon2opensim/model/lower_limb.osim"
    file_model_moco                 = "/home/lisca/biomec/sim/src/opensim-core/OpenSim/Examples/Moco/example3DWalking/subject_walk_armless.osim"
    file_model_vicon2opensim2moco   = "/home/lisca/biomec/sim/src/athlete/matlab/vicon2opensim/model/vicon2opensim2moco.osim"

    converter = Vicon2OpenSim2MocoModelConverter(
        file_model_opensim_matlab_tools,
        file_model_moco,
        file_model_vicon2opensim2moco,)

    converter.adjust_marker_used()
    converter.adjust_marker_positions()
    converter.adjust_joints()
    converter.save_to_file()
