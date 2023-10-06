import mill._
import mill.scalalib._
import mill.scalalib.publish._
import java.util.jar.Attributes.Name.{IMPLEMENTATION_VERSION => ImplementationVersion}

import coursier.maven.MavenRepository

import scala.sys.process.Process
import scala.util.{Failure, Success, Try}

/** Base trait for build modules. */
trait CommonModule extends SbtModule {
  def deployLocal(assembly: PathRef, jarName:String) = {
    os.makeDir.all(os.pwd / Symbol("jars"))
    println(s"Copying artifact ${assembly.path} to jars / $jarName")
    os.copy(assembly.path, os.pwd / Symbol("jars") / jarName, replaceExisting = true)
  }

  override def repositories: Seq[coursier.Repository] = super.repositories ++ Seq(
    MavenRepository("https://oss.sonatype.org/content/repositories/public"),
    MavenRepository("https://oss.sonatype.org/content/repositories/snapshots"),
    MavenRepository("https://jcenter.bintray.com/"),
    MavenRepository("https://broadinstitute.jfrog.io/artifactory/libs-snapshot/")
  )
}

/** A base trait for versioning modules. */
trait ReleaseModule extends JavaModule {
  /** Execute Git arguments and return the standard output. */
  private def git(args: String*): String = os.proc("git", args).call().out.text().trim

  /** Get the commit hash at the HEAD of this branch. */
  private def gitHead: String = git("rev-parse", "HEAD")

  /** Get the commit shorthash at the HEAD of this branch .*/
  private def shortHash: String = gitHead.take(7)

  /** The current tag of the currently checked out commit, if any. */
  def currentTag: Try[String] = Try(git("describe", "--exact-match", "--tags", "--always", gitHead))

  /** The hash of the last tagged commit. */
  private def hashOfLastTag: Try[String] = Try(git("rev-list", "--tags", "--max-count=1"))

  /** The last tag of the currently checked out branch, if any. */
  def lastTag: Try[String] = hashOfLastTag match {
    case Success(hash) => Try(git("describe", "--abbrev=0", "--tags", "--always", hash))
    case Failure(e)    => Failure(e)
  }

  /** If the Git repository is left in a dirty state. */
  private def dirty: Boolean = git("status", "--porcelain").nonEmpty

  /** The implementation version. */
  private def implementationVersion = T.input {
    val prefix: String = (currentTag, lastTag) match {
      case (Success(_currentTag), _)       => _currentTag
      case (Failure(_), Success(_lastTag)) => _lastTag + "-" + shortHash
      case (_, _)                          => shortHash
    }
    prefix + (if (dirty) "-dirty" else "")
  }

  /** The version string `Target`. */
  def version = T input { println(implementationVersion()) }

  /** The JAR manifest. */
  override def manifest = T { super.manifest().add(ImplementationVersion.toString -> implementationVersion()) }
}


object tools extends CommonModule with PublishModule with ReleaseModule {
  def scalaVersion = "2.13.8"
  override def millSourcePath = super.millSourcePath / os.up
  override def mainClass = Some("com.fulcrumgenomics.sv.cmdline.SvMain")
  override def artifactName = "fgsv"
  def gitHash = Process("git rev-parse --short HEAD").lineStream.head
  def publishVersion = s"0.1.0-${gitHash}-SNAPSHOT"
  def pomSettings = PomSettings(
    description = artifactName(),
    organization = "com.fulcrumgenomics",
    url = "https://github.com/fulcrumgenomics/fgsv",
    licenses = Seq(License.MIT),
    versionControl = VersionControl.github("fulcrumgenomics","fgsv"),
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

  override def ivyDeps = Agg(
    ivy"org.scala-lang:scala-compiler:${scalaVersion()}",
    ivy"com.fulcrumgenomics:fgbio_2.13:2.1.0".excludeOrg(orgsToExclude:_*)
  )

  object test extends Tests {
    override def ivyDeps = Agg(ivy"org.scalatest::scalatest:3.1.0")
    override def testFramework = "org.scalatest.tools.Framework"

    // run mill tools.test.singleTest com.fulcrumgenomics.sv.x.y.x.TestClassName
    def singleTest(args: String*) = T.command {
      super.runMain("org.scalatest.run", args: _*)
    }
  }

  override def javacOptions = Seq("-source", "1.8", "-target", "1.8")

  def deployLocal = T { super.deployLocal(assembly(), "fgsv.jar")  }
}
