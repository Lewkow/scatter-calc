
package scatter

import potential.Potential

class Scatter(energyC: Double, 
              projectileC: String, 
              targetC: String,
              potential_typeC: String) {

  import math.sqrt
  val energy: Double = energyC
  val coll = new Collision(projectileC, targetC)
  val projectile: String = projectileC
  val target: String = targetC
  val potential_type: String = potential_typeC
  val wavenumber: Double = math.sqrt(2.0d*coll.reduced_mass*energy)
  val phase = new Phaseshift(energy, potential_type, coll)

  def printValues() {
    val grid_start = phase.grid.start_pos 
    val grid_end = phase.grid.end_pos 
    println(s"Energy:       $energy [eV]")
    println(s"Grid Start:   $grid_start [a0]")
    println(s"Grid End:     $grid_end [a0]")
  }

}

class Collision(projectile: String,
                target: String) {

  val masses = Map("H"->1.00782503207d,
                   "He3"->3.0160293191d,
                   "He4"->4.00260325415d,
                   "O"->15.99491461956d)
  val projectile_mass: Double = masses(projectile)
  val target_mass: Double = masses(target)
  val reduced_mass: Double = reducedMass(projectile_mass,target_mass)
  printValues()

  def reducedMass(proj: Double, targ: Double) = {
    val denom: Double = proj + targ
    val reduced_mass: Double = proj*targ/denom
    reduced_mass
  }

  def printValues() {
    println(s"Projectile:   $projectile") 
    println(s"Mass:         $projectile_mass [u]")
    println(s"Target:       $target")
    println(s"Mass:         $target_mass [u]")
    println(s"Reduced Mass: $reduced_mass [u]")
  }

}

class Phaseshift(energy: Double, potential_type: String, coll: Collision) {
  val MAX_GRID_SIZE = 100000
  val pot = new Potential(potential_type)
  val grid = new Grid(energy, pot, coll)
}

class Grid(energyC: Double, pot: Potential, coll: Collision) {
  val energy: Double = energyC
  val start_pos: Double = grid_start_finder() 
  val end_pos: Double = grid_end_finder()
  val dr_grid: Double = grid_dr_finder()
  val N_grid: Int = grid_N_finder() 
  val pos_grid: Array[Double] = build_pos_grid()
  val pot_grid: Array[Double] = build_pot_grid()

  def grid_dr_finder(): Double = {
    val k_max: Double = math.sqrt(2.0d*coll.reduced_mass*(energy+pot.depth))
    val h_max: Double = (1.0d/k_max)/20.0d
    var step: Double = 7.0e-4
    if (7.0e-4 > h_max) {
      step = h_max
    }
    step
  }  

  def grid_N_finder() = {
    val N_grid: Int = ((end_pos - start_pos)/dr_grid).toInt + 1
    N_grid
  }

  def grid_end_finder() = {
    val lim: Double = 1e-8
    val dr: Double = 0.1d
    var r: Double = start_pos
    var pot1: Double = pot.get_potential(1.0d)
    var pot2: Double = pot1
    var dpot: Double = 0.0d
    var going : Boolean = true
    var i: Int = 1
    var r_end: Double = start_pos

    while (going == true) {
      r = start_pos + i*dr
      pot2 = pot.get_potential(r)
      dpot = math.abs(pot2-pot1)
      if (dpot < lim) {
        r_end = r
        going = false
      }
      pot1 = pot2
      i = i + 1
    }
    r_end
  }

  def grid_start_finder() = {
    var i: Int = 0
    var j: Int = 0 
    val start: Double = 10.0d
    val dr: Double = 0.0001
    val lim: Double = energy*1.3d
    var going: Boolean = true
    var r: Double = start

    while (going == true) {
      r = start - j.toDouble*dr
      var pot_r: Double = pot.get_potential(r)
      if ( (pot_r >= lim) | (r <= 1.0e-5) ) {
        going = false 
      } else {
        j = j + 1
      }
    }
    r
  }

  def build_pos_grid() = {
    val z = new Array[Double](N_grid)    
    var i: Int = 0
    for (i <- 0 to z.length-1) {
      z(i) = start_pos + i.toDouble*dr_grid
    }
    z 
  }

  def build_pot_grid() = {
    val z = new Array[Double](N_grid)
    var i: Int = 0
    for (i <- 0 to z.length-1) {
      var x: Double = start_pos + i.toDouble*dr_grid
      z(i) = pot.get_potential(x)
    }   
    z
  }

}

