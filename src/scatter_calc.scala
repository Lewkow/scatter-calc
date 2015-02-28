package runtime

import scatter.Scatter

object Main extends App{
  val test_scatter = new scatter.Scatter(10.0d, "H", "H")
  test_scatter.printValues()
}

// object scatter_calc {

//   def main(args: Array[String]) {
//     val test_scatter = Scatter(10.0d, "H", "H")
//     test_scatter.printValues()
//   }

// }

