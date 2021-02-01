package cn.edu.pku.asic.storage.dggs.h3;

import com.uber.h3core.H3Core;
import com.uber.h3core.exceptions.DistanceUndefinedException;
import com.uber.h3core.exceptions.LineUndefinedException;
import com.uber.h3core.exceptions.PentagonEncounteredException;
import com.uber.h3core.util.GeoCoord;

import java.io.IOException;
import java.util.List;

/**
 * =============================================================================================
 *
 * @ProjectName dggs
 * @Package cn.edu.pku.asic.storage.dggs.h3
 * @ClassName test
 * @Version 1.0.0
 * @Author Guoliang PU
 * @Date 2021/2/1 19:12
 * @See =============================================================================================
 * @Description
 */
public class test {

    public static void main(String[] args) throws IOException, DistanceUndefinedException, PentagonEncounteredException, LineUndefinedException {
        //将经纬度编码为六边形地址
        H3Core h3 = H3Core.newInstance();
        double lat = 31.920523;
        double lng = 119.168182;
        //H3分辨率
        int res = 6;
        String hexAddr = h3.geoToH3Address(lat, lng, res);
        System.out.println("H3地址：" + hexAddr);

        //H3地址转换为经纬度
        List<GeoCoord> geoCoords = h3.h3ToGeoBoundary(hexAddr);
        System.out.println("六边形顶点经纬度：" + geoCoords);

        //geoToH3：经纬度转换成H3索引
        System.out.println("H3索引Long类型：" + h3.geoToH3(lat, lng, 6));

        //查找索引的质心
        System.out.println("H3索引质心,经纬度坐标类型：" + h3.h3ToGeo(617847554087583743L));

        //查找索引的边界
        System.out.println("六边形顶点经纬度：" + h3.h3ToGeoBoundary(604336920381620223L));

        //根据提供的起点和终点返回单向边缘H3索引
        //起点镇江市
        double originLat = 31.977022;
        double originLng = 119.160861;

        //起终点经纬度转换为H3索引
        Long originH3Index = h3.geoToH3(originLat, originLng, 9);
        Long destH3Index = h3.geoToH3(lat, lng, 9);

        //查找相邻索引,将k=1定义为所有相邻索引
        System.out.println("origin相邻索引：" + h3.kRing(originH3Index, 1));

        //两索引是否为邻居
        boolean isNeighbors = h3.h3IndexesAreNeighbors(originH3Index, 617847554084700159L);
        System.out.println("两索引是否为邻居：" + isNeighbors);

        //计算单向边缘H3索引,单向边缘允许对从一个小区到相邻小区的有向边缘进行编码
        System.out.println("单向边缘H3索引:" + h3.getH3UnidirectionalEdge(originH3Index, 617847554084700159L));

        //返回两个索引之间的网格单元格距离
        System.out.println("两个索引之间的网格单元格距离：" + h3.h3Distance(originH3Index, destH3Index));

        //返回空心六角环
        System.out.println("周围相邻空心六角环：" + h3.hexRing(originH3Index, 1));

        //返回两个索引之间的索引行
        System.out.println("两个索引之间的索引行：" + h3.h3Line(originH3Index, destH3Index));
    }
}