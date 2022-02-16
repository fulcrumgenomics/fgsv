import java.text.SimpleDateFormat
import java.util.Date

import mill._
import mill.scalalib._
import mill.scalalib.publish.{Developer, License, PomSettings, SCM}
import java.util.jar.Attributes.Name.{IMPLEMENTATION_VERSION => ImplementationVersion}

import ammonite.ops._
import coursier.maven.MavenRepository

import scala.sys.process.Process
import scala.util.{Failure, Success, Try}

/** Base trait for build modules. */
trait CommonModule extends SbtModule {
  def deployLocal(assembly: PathRef, jarName:String) = {
    mkdir(pwd / 'jars)
    println(s"Copying artifact ${assembly.path} to jars / $jarName")
    cp.over(assembly.path, pwd / 'jars / jarName)
  }

  override def repositories: Seq[coursier.Repository] = super.repositories ++ Seq(
    MavenRepository("https://oss.sonatype.org/content/repositories/public"),
    MavenRepository("https://oss.sonatype.org/content/repositories/snapshots"),
    MavenRepository("https://jcenter.bintray.com/"),
    MavenRepository("https://artifactory.broadinstitute.org/artifactory/libs-snapshot/")
  )
}

/** A base trait for versioning modules. */
trait ReleaseModule extends JavaModule {
  /** Execute Git arguments and return the standard output. */
  private def git(args: String*): String = %%("git", args)(pwd).out.string.trim

  /** Get the commit hash at the HEAD of this branch. */
  private def gitHead: String = git("rev-parse", "HEAD")

  /** Get the commit shorthash at the HEAD of this branch .*/
  private def shortHash: String = gitHead.take(7)

  /** If the Git repository is left in a dirty state. */
  private def dirty: Boolean = git("status", "--porcelain").nonEmpty

  private def today: String = new SimpleDateFormat("yyyyMMdd").format(new Date())

  /** The implementation version. */
  private def implementationVersion = T input {
    val prefix = s"${today}-${shortHash}"
    if (dirty) s"${prefix}-dirty" else prefix
  }

  /** The version string `Target`. */
  def version = T input { println(implementationVersion()) }

  /** The JAR manifest. */
  override def manifest = T { super.manifest().add(ImplementationVersion.toString -> implementationVersion()) }
}


object tools extends CommonModule with PublishModule with ReleaseModule {
  def scalaVersion = "2.13.3"
  def millSourcePath = super.millSourcePath / ammonite.ops.up
  def mainClass = Some("com.fulcrumgenomics.sv.cmdline.SvMain")
  def artifactName = "fgsv"
  def gitHash = Process("git rev-parse --short HEAD").lineStream.head
  def publishVersion = s"0.0.1-${gitHash}-SNAPSHOT"
  def pomSettings = PomSettings(
    description = artifactName(),
    organization = "com.fulcrumgenomics",
    url = "https://github.com/fulcrumgenomics/fgsv",
    licenses = Seq(License("MIT License", "http://www.fulcrumgenomics.com")),
    scm = SCM(
      "git://github.com:fulcrumgenomics/fgsv.git",
      "scm:git://github.com:fulcrumgenomics/fgsv.git"
    ),
    developers = Seq(
      Developer("tfenne", "Tim Fennell", "https://github.com/tfenne"),
      Developer("nh13", "Nils Homer", "https://github.com/nh13"),
    )
  )

  private val orgsToExclude = Seq(
    "org.apache.ant",
    "gov.nih.nlm.ncbi",
    "org.testng",
    "com.google.cloud.genomics"
  )

  def ivyDeps = Agg(
    ivy"org.scala-lang:scala-compiler:${scalaVersion()}",
    ivy"com.fulcrumgenomics:fgbio_2.13:1.5.0".excludeOrg(orgsToExclude:_*),
	ivy"org.xerial.snappy:snappy-java:1.1.8.4"
    )

  object test extends Tests {
    def ivyDeps = Agg(ivy"org.scalatest::scalatest:3.1.0")
    def testFramework = "org.scalatest.tools.Framework"

    // run mill tools.test.singleTest com.fulcrumgenomics.sv.x.y.x.TestClassName
    def singleTest(args: String*) = T.command {
      super.runMain("org.scalatest.run", args: _*)
    }
  }

  def javacOptions = Seq("-source", "1.8", "-target", "1.8")

  def deployLocal = T { super.deployLocal(assembly(), "fgsv.jar")  }
}
