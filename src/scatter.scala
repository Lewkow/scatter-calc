package scatter

class Scatter(energyC: Double, 
              projectileC: String, 
              targetC: String) {

  val masses = Map("H"->1.00782503207d,
                   "He3"->3.0160293191d,
                   "He4"->4.00260325415d,
                   "O"->15.99491461956d)
  import math.sqrt
  val energy: Double = energyC
  val projectile: String = projectileC
  val target: String = targetC
  val projectile_mass: Double = masses(projectile)
  val target_mass: Double = masses(target)
  val reduced_mass: Double = reducedMass(projectile_mass,target_mass)
  val wavenumber: Double = math.sqrt(2.0d*reduced_mass*energy)

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

  def setup_phaseshifts() {
    var phase = new Phaseshift()
  }

}

class Phaseshift() {
  val MAX_GRID_SIZE = 100000
  var grid_step_size: Double = 0.0d
  var grid_start: Double = 0.0d
  var grid_end: Double = 0.0d

  import potential.Potential
  val tester = new Potential("LJ")
  var i = 0
  for (i <- 0 to 30) {
    var r: Double = i.toDouble + 0.001d
    var v = tester.get_potential(r)
    println(s"V($r) = $v") 
  }

  def grid_start_finder(energy: Double) = {
    0.0d
  }

}

