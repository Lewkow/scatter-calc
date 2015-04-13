package runtime

import java.io._
import scatter.Scatter

object scatter_calc {

  def main(args: Array[String]) {
    // val test_scatter = new Scatter(100.0d, "H", "H", "LJ")
    // test_scatter.printValues()
    test_differential_cross_section_file
  }

  def get_header(calc_type: String, 
                 energy: Double,
                 proj: String,
                 targ: String,
                 pot_type: String):String = {
    val a = "-".toString
    calc_type + a + energy.toString + a + proj + a + targ + a + pot_type
  }

  def test_differential_cross_section_file() {
    val energy: Double = 10.0d
    val proj: String = "H"
    val targ: String = "H"
    val pot: String  = "LJ"
    val test_scatter = new Scatter(energy, proj, targ, pot)
    val N: Int       = 100
    val dt: Double   = math.Pi/(N-1).toDouble
    val theta        = new Array[Double](N)
    val dcs          = new Array[Double](N)
    var i: Int       = 0
    for (i <- 0 to N-1) {
      theta(i) = i*dt
      dcs(i)   = test_scatter.cross_sections.differential_cross_section(theta(i))
    }
    val header:String = get_header("dcs",energy,proj,targ,pot)
    print_xy_to_file("./data/scatter_calc_out.dat", header, theta, dcs)
  }

  def print_xy_to_file(filename: String, header: String, x: Array[Double], y: Array[Double]) {
    if (x.length != y.length) {
      println("X and Y lengths NOT equal!")
    } else {
      val fw = new FileWriter(filename, true)
      fw.write(header+"\n")
      var i: Int = 0 
      try {
        for (i <- 0 to x.length-1) {
          val s: String = List(x(i), y(i)) mkString(", ")
          val n: String = "\n"
          val z: String = s + n
          fw.write(z)
        }
      }
      finally fw.close()
    }
  }



}

