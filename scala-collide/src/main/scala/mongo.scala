package mongo

import com.mongodb.casbah.Imports._

class Mongo() {
  val mongoClient = MongoClient("localhost", 27017)
  val db = mongoClient("scatter")
  
}