
package spark_scatter

// Imported classes/packages
import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._
import org.apache.spark.SparkConf
import org.apache.spark.sql._
import org.apache.spark.rdd._

import scala.collection.mutable.ArrayBuffer


class SparkScatter {

  val conf     = new SparkConf().setAppName("CIS Views")
  val sc       = new SparkContext(conf)

  val x = 1 to 100 toArray
  val y = easy_add(sc, x)
  println("Sum: " + y.toString) 

  def easy_add(sc: SparkContext, x: Array[Int]): Int = {
    val r = sc.parallelize(x)
    val y = r.reduce(_+_)
    y
  }

}
