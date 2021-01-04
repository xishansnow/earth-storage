package cn.edu.pku.coe.asic.earthstorage.indexing;

import org.apache.hadoop.util.IndexedSortable;

/**
 *=============================================================================================
 * @ProjectName    EarthStorage 
 * @Package        cn.edu.pku.coe.asic.earthStorage.indexing
 * @ClassName      GridSortable
 * @Version        1.0.0
 * @Author         Guoliang PU
 * @Date           2021/1/2 20:46
 * @See
 *=============================================================================================
 * @Description    A class for comparing and sorting points based on the grid coding schema.
 */
public class GridSortable  implements IndexedSortable {


    @Override
    public void swap(int i, int i1) {

    }

    @Override
    public int compare(int i, int i1) {
        return 0;
    }
}
