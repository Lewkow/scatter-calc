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

}