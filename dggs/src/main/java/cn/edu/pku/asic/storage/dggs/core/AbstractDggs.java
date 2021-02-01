package cn.edu.pku.asic.storage.dggs.core;

import cn.edu.pku.asic.storage.dggs.sphere.SpherePoint;

import java.io.Externalizable;
import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/**
 *=============================================================================================
 * @ProjectName    EarthStorage 
 * @Package        cn.edu.pku.asic.storage.dggs
 * @ClassName      AbstractDggs
 * @Version        1.0.0
 * @Author         Guoliang PU
 * @Date           2021/1/29 13:41
 * @See
 *=============================================================================================
 * @Description    Parent of all DGG System.            
 */
public abstract class AbstractDggs<TID extends AbstractCellId,TCELL extends AbstractCell,TCELLUNION extends AbstractCellUnion> implements Externalizable {
    /**
     *  Metadata annotates the subclass of AbstractDGGS, to describe it's basic feature descriptions
     */
    @Target(ElementType.TYPE)
    @Retention(RetentionPolicy.RUNTIME)
    @interface Metadata
    {
        // DGGS基本属性
        String dggsType() default "NONE";              // type of DGGS,
        // STANDARD,CUSTOM,SUPERFUND,PLANETRISK etc.
        // STANDARD expression == ProjectionName + Aperture + Topology（Abbreviated）
        // Examples: ISEA4T means ISEA projection, Pure Aperture of 4, and Triangle topology.
        // FULLER7H means FULLER projection, Pure Aperture of 7, and Hexagon topology
        String dggsTopology()default "NONE";          // value in {HEXAGON,TRIANGLE,DIAMOND,SQUARE}
        String dggsProj()default "NONE";              // value in {ISEA,FULLER,MERCATOR,QTM,EARPIH,SQT}

        //剖分层级与剖分方案
        int dggsResSpec()default 0;              // specified DGG resolution levels, value in [0,35]
        // if dggs_type is SUPERFUND  then 0 ≤ v ≤ 9;
        // if dggs_aperture_type is SEQUENCE then 0 ≤ v ≤ n, where n is  the length of  dggs_aperture_sequence
        // DGGS的分辨率层级数量，取值范围：0到35级
        String dggsResSpecifyType() default "NONE";    // how is the DGGS resolution specified? value in {SPECIFIED,CELL_AREA,INTERCELL_DISTANCE}
        // 网格分辨率的设置方法，{按指定大小设置分辨率，按cell面积设置分辨率，按cell间距离设置分辨率}
        double dggsResSpecifyArea() default 100;    // disired cell area,value in [1.0,..)
        // 如果按照CELL_AREA方案，设置分辨率值，最小值为1.0
        double dggsResSpecifyIntercellDistance()default 100;
        // desired intercell distance (measured on the plane),value in [1.0,..)
        // 如果按照INTERCELL_DISTANCE方案，则设置距离值，最小值为1.0
        /*     孔径类型和参数   */
        String dggsApertureType()default "PURE";      // value in {PURE，MIXED43，SEQUENCE}
        int dggsAperture()default 4;             // aperture of DGGS if type is PURE, integer of 3,4,7 or 9
        String dggsApertureSequence()default "44444444444444444444444444444444444";  // the DGGS aperture sequence if aperture type is SEQUENCE
        // string of 3’s, 4’s,and 7’s in any order, like "333333333333"
        int dggsNumAperture4Res()default 9;      // number of aperture 4 resolutions in a mixed aperture sequence, value in [0,35]

        /*    参考椭球体       */
        String proj_datum()default "NONE";            // desired earth radius datum, value in {WGS84_AUTHALIC_SPHERE,WGS84_MEAN_SPHERE, CUSTOM_SPHERE}

        /*  剖分原点的方案     */
        String dggsOrientSpecifyType()default "NONE"; // how is the DGG orientation specified? value in {RANDOM,SPECIFIED,REGION_CENTER}
        // if SPECIFIED:
        double dggsVert0Azimut()default 0.0;       // azimuth from icosahedron  vertex 0 to vertex 1 (degrees)，value in [-90.0 ≤ v ≤ 90.0]
        double dggsVert0Lat()default 0.0;          // latitude from icosahedron  vertex 0 to vertex 1 (degrees)，value in [-90.0 ≤ v ≤ 90.0]
        double dggsVert0Long()default 0.0;         // longitude  from icosahedron  vertex 0 to vertex 1 (degrees)
        // if REGION_CENTER：
        double regionCenterLat()default 0.0;       // latitude of study region (degrees)
        double regionCenterLong()default 0.0;      // longitude of study region (degrees)

        String description() default "NONE";           // free text description
    }



    public SpherePoint origin(){            // return origin point's coordinate of DGG system

        SpherePoint origin = null;

        return(origin);
//        Metadata mt = getClass().getAnnotation(Metadata.class);
//        Point origin;
//        origin.mt.regionCenterLat(),mt.regionCenterLong()
    }

    public String projDatum(){              // return dataum of DGG system
        return "NONE";
    }
    public int resolutionSpec(){            // return DGG system's resolution levels, value in [0,35]
        return 4;
    }
    public String topologyType(){           // return
        return "NONE";
    }
    public String resolutionSpecifyType(){  // how is the DGGS resolution specified? value in {SPECIFIED,CELL_AREA,INTERCELL_DISTANCE}
        return "NONE";
    }
    public String apertureType(){           // value in {PURE，MIXED43，SEQUENCE}
        return "NONE";
    }
    public String apertureSequence(){             // aperture sequence in string format
        return "NONE";
    }
}
