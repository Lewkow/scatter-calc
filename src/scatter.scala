
package scatter

class Scatter(energyC: Double, 
              projectileC: String, 
              targetC: String,
              potential_typeC: String) {

  val masses = Map("H"->1.00782503207d,
                   "He3"->3.0160293191d,
                   "He4"->4.00260325415d,
                   "O"->15.99491461956d)
  import math.sqrt
  val energy: Double = energyC
  val projectile: String = projectileC
  val target: String = targetC
  val potential_type: String = potential_typeC
  val projectile_mass: Double = masses(projectile)
  val target_mass: Double = masses(target)
  val reduced_mass: Double = reducedMass(projectile_mass,target_mass)
  val wavenumber: Double = math.sqrt(2.0d*reduced_mass*energy)
  val phase = new Phaseshift(potential_type)



  def reducedMass(proj: Double, targ: Double) = {
    val denom: Double = proj + targ
    val reduced_mass: Double = proj*targ/denom
    reduced_mass
  }

  def printValues() {
    println(s"Projectile: $projectile Mass: $projectile_mass [u]")
    println(s"Target: $target Mass: $target_mass [u]")
    println(s"Energy: $energy [eV]")
    println(s"Reduced Mass: $reduced_mass [u]")
  }

}

class Phaseshift(potential_type: String) {
  val MAX_GRID_SIZE = 100000
  // var grid_step_size: Double = 0.0d
  var grid_start: Double = 0.001d
  var grid_end: Double = 10.0d
  val grid = new Grid(grid_start, 
                  grid_end, 
                  10,
                  // MAX_GRID_SIZE,
                  potential_type)
  // val pos_grid = grid.build_pos_grid()
  // val pot_grid = grid.build_pot_grid(potential_type)

  def grid_start_finder(energy: Double) = {
    0.0d
  }

}

class Grid(startC: Double, 
           endC: Double, 
           N_gridC: Int,
           potential_typeC: String) {
  val start_pos: Double = startC
  val end_pos: Double = endC
  val N_grid: Int = N_gridC
  val dx: Double = (end_pos - start_pos)/(N_grid.toDouble - 1.0d)
  val potential_type: String = potential_typeC
  val pos_grid = build_pos_grid()
  val pot_grid = build_pot_grid()

  def build_pos_grid() = {
    val z = new Array[Double](N_grid)    
    var i: Int = 0
    for (i <- 0 to z.length-1) {
      z(i) = start_pos + i.toDouble*dx
    }
    z 
  }

  def build_pot_grid() = {
    import potential.Potential
    val pot = new Potential(potential_type)
    val z = new Array[Double](N_grid)
    var i: Int = 0
    for (i <- 0 to z.length-1) {
      var x: Double = start_pos + i.toDouble*dx
      z(i) = pot.get_potential(x)
    }   
    z
  }



}

