package potential

class Potential(potential_methodC: String) {
  val potential_method: String = potential_methodC
  val depth = get_depth(true)
  val depth_pos = get_depth(false)
  // print_values()

  // pass interparticle distance [u]
  def get_potential(r: Double): Double = {
    if (potential_method == "LJ") {
      lennard_jones(r, 1.0d, 1.0d)
    }
    else {
      1.0d
    }
  }

  def get_depth(depth_mode: Boolean): Double = {
    val r_start: Double = 0.5d
    val r_end: Double = 20.0d
    val dr: Double = 0.001 
    var r: Double = r_start
    var i: Int = 0
    var pot: Double = 0
    var deepest: Double = 0.0d
    var deep_pos: Double = 0.0d

    while (r < r_end) {
      r = r_start + i.toDouble*r
      pot = get_potential(r)
      if (pot < deepest) {
        deepest = pot
        deep_pos = r
      }
      i = i + 1
    }
    if (depth_mode == true) {
      deepest
    }
    else {
      deep_pos
    }
  }

  def lennard_jones(r: Double, 
                    well_depth: Double,
                    wall_pos: Double): Double = {
    if (r <= 0.0d) {
      0.0d
    }
    else {
      val x = wall_pos/r
      4.0d*well_depth*(math.pow(x,12.0)-math.pow(x,6.0))
    }
  }

  def print_values() = {
    println(f"Depth:        $depth%1.3f [eV]")
    println(f"Depth Pos:    $depth_pos%1.3f [a0]")
  }
}
