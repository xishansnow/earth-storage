package cn.edu.pku.coe.asic.earthstorage;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/**
 *=============================================================================================
 * @ProjectName    EarthStorage 
 * @Package        cn.edu.pku.coe.asic.earthStorage
 * @ClassName      JavaTest
 * @Version        1.0.0
 * @Author         Guoliang PU
 * @Date           2021/1/3 9:24
 * @See
 *=============================================================================================
 * @Description    A main class for testing Java grammar            
 */

//Annotation：Java中利用反射机制形成的注释机制
//Annotation: 不影响程序运行，但会对编译器等辅助工具产生影响
//Annotation本身是一个接口，@interface不是接口，而是Annotation接口的一个实现
//Annotation的定义方法 public @interface XXXX{}
//Annotation可以定义自己的成员变量，也可以不含任何成员变量（如上行所示），只用于编程时的标记和提醒
//用@Target来设置所定义Annotation适用的程序元素类型（如：类型、构造函数、成员变量、方法、参数、包等
//用@Retention来设置所定义Annotation的有效范围（如：仅在源码中可使用，在编译的类文件中可使用，在运行时也可使用）
//Annotation的使用方法，在待修饰的程序元素前，输入@XXXX来引用
//典型例子：
//  变量的注释：  @XXXX String m_strMember;
//  方法的注释：  @XXXX public Print(int id, String name)
//  参数的注释    public Print(@XXXX int id, @XXXX String name)

//Annotation的访问：当@Retention设置为RUNTIME时，可以在程序中访问程序元素的Annotation信息（具体待用时在细致了解）


@Target(ElementType.FIELD)
@Retention(RetentionPolicy.SOURCE)
@interface  Preferred {}

@Target(ElementType.PARAMETER)
@Retention(RetentionPolicy.RUNTIME)
@interface Required {}

//class used to test reflection
class TextField{
    int a,b;

    TextField() {

    }


    Print(){

    }

}



public class JavaTest {
    public static void main(String[] args) {
        System.out.println("Hello，I'm EarthStorage!");
        TextField tf = new TextField();
        Class textFieldC = tf.getClass();
        System.out.println(textFieldC);
    }

}


