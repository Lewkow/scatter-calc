package potential

class Potential(potential_methodC: String) {
  val potential_method: String = potential_methodC

  // pass interparticle distance [u]
  def get_potential(r: Double): Double = {
    if (potential_method == "LJ") {
      lennard_jones(r, 1.0d, 1.0d)
    }
    else {
      1.0d
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
}
