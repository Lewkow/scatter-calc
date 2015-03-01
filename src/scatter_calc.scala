package runtime

import scatter.Scatter

object scatter_calc {

  def main(args: Array[String]) {
    val test_scatter = new Scatter(10.0d, "H", "H")
    test_scatter.printValues()
    test_scatter.setup_phaseshifts()
  }

}

