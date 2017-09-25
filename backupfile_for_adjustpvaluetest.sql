-- MySQL dump 10.13  Distrib 5.7.19, for osx10.10 (x86_64)
--
-- Host: localhost    Database: adjustpvaluetest
-- ------------------------------------------------------
-- Server version	5.7.19

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `density`
--

DROP TABLE IF EXISTS `density`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `density` (
  `snp_id` varchar(45) NOT NULL,
  `test_id` varchar(45) DEFAULT NULL,
  `lower` int(11) DEFAULT NULL,
  `upper` int(11) DEFAULT NULL,
  `naive_pv` double DEFAULT NULL,
  `adjusted_pv` double DEFAULT NULL,
  PRIMARY KEY (`snp_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `density`
--

LOCK TABLES `density` WRITE;
/*!40000 ALTER TABLE `density` DISABLE KEYS */;
INSERT INTO `density` VALUES ('rs11','cd',2,5,0.00000008,0.6),('rs12','r',6,8,0.00000007,0.8),('rs13','d',10,11,0.000000009,0.6),('rs14','a',12,14,0.0000000018,0.2),('rs15','d',16,19,0.00000007,0.2);
/*!40000 ALTER TABLE `density` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `genotype2`
--

DROP TABLE IF EXISTS `genotype2`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `genotype2` (
  `snp_id` varchar(45) NOT NULL,
  `snp_genotype` text,
  PRIMARY KEY (`snp_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `genotype2`
--

LOCK TABLES `genotype2` WRITE;
/*!40000 ALTER TABLE `genotype2` DISABLE KEYS */;
INSERT INTO `genotype2` VALUES ('a','021102'),('b','101021'),('c','000110'),('d','210021'),('e','100221'),('f','001110'),('g','221012'),('h','010221'),('i','201012'),('j','020111'),('k','030112'),('l','102211'),('m','021100'),('n','021021'),('o','011211'),('p','102012'),('q','200021'),('r','303112'),('s','012021'),('t','013200');
/*!40000 ALTER TABLE `genotype2` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `phenotype`
--

DROP TABLE IF EXISTS `phenotype`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `phenotype` (
  `id` int(11) NOT NULL,
  `phenotype_string` text,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `phenotype`
--

LOCK TABLES `phenotype` WRITE;
/*!40000 ALTER TABLE `phenotype` DISABLE KEYS */;
INSERT INTO `phenotype` VALUES (1,'shamrock,dandelion,magenta,denim,shamrock,magenta');
/*!40000 ALTER TABLE `phenotype` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Table structure for table `phenotype2`
--

DROP TABLE IF EXISTS `phenotype2`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `phenotype2` (
  `id` int(11) NOT NULL,
  `phenotype_string` text,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `phenotype2`
--

LOCK TABLES `phenotype2` WRITE;
/*!40000 ALTER TABLE `phenotype2` DISABLE KEYS */;
INSERT INTO `phenotype2` VALUES (1,'130210');
/*!40000 ALTER TABLE `phenotype2` ENABLE KEYS */;
UNLOCK TABLES;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2017-09-25 12:57:42
