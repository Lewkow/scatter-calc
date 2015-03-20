package runtime

import scatter.Scatter

object scatter_calc {

  def main(args: Array[String]) {
    val test_scatter = new Scatter(100.0d, "H", "H", "LJ")
    test_scatter.printValues()
  }

}

