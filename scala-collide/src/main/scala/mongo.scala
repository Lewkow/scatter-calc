package mongo

import com.mongodb.casbah.Imports._

class Mongo() {
  def truncateAt(n: Double, p: Int): Double = { 
    val s = math pow (10, p); (math floor n * s) / s 
  }

  val float_precision: Int = 5

  val mongoClient = MongoClient("localhost", 27017)
  val db = mongoClient("scatter")
  
  def write_tcs(collision_id: String, energy: Double, tcs: Double) {
    val coll = db("tcs")
    val id = "tcs_" + collision_id + "_" + truncateAt(energy, float_precision).toString
    val prev_doc = MongoDBObject("_id" -> id)
    val prev = coll.findOne(prev_doc)
    if (prev.isEmpty == true) {
      val doc = MongoDBObject("_id" -> id, 
                              "energy" -> truncateAt(energy, float_precision),
                              "tcs" -> truncateAt(tcs, float_precision))
      coll.insert(doc)
    }
  }

  def write_dcs(collision_id: String, energy: Double, angle: Double, dcs: Double) {
    val coll = db("dcs")
    val id = "dcs_" + collision_id + "_" + truncateAt(energy, float_precision).toString +
             "_" + truncateAt(angle, float_precision).toString
    val prev_doc = MongoDBObject("_id" -> id)
    val prev = coll.findOne(prev_doc)
    if (prev.isEmpty == true) {
      val doc = MongoDBObject("_id" -> id, 
                              "energy" -> truncateAt(energy, float_precision),
                              "angle" -> truncateAt(angle, float_precision),
                              "dcs" -> truncateAt(dcs, float_precision))
      coll.insert(doc)
    }
  }
}

