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

  override def repositoriesTask: Task[Seq[coursier.Repository]] = T.task {
    super.repositoriesTask() ++ Seq(
      MavenRepository("https://oss.sonatype.org/content/repositories/public"),
      MavenRepository("https://oss.sonatype.org/content/repositories/snapshots"),
      MavenRepository("https://jcenter.bintray.com/"),
      MavenRepository("https://broadinstitute.jfrog.io/artifactory/libs-snapshot/")
    )
  }

  /** All Scala compiler options for this package. */
  override def scalacOptions: T[Seq[String]] = T {
    Seq(
      "-opt:inline:com.fulcrumgenomics.**", // Turn on the inliner.
      "-opt-inline-from:com.fulcrumgenomics.**", // Tells the inliner that it is allowed to inline things from these classes.
      "-Yopt-log-inline", "_", // Optional, logs the inliner activity so you know it is doing something.
      "-Yopt-inline-heuristics:at-inline-annotated", // Tells the inliner to use your `@inliner` tags.
      "-opt-warnings:at-inline-failed", // Tells you if methods marked with `@inline` cannot be inlined, so you can remove the tag.
      // The following are sourced from https://nathankleyn.com/2019/05/13/recommended-scalac-flags-for-2-13/
      "-deprecation", // Emit warning and location for usages of deprecated APIs.
      "-explaintypes", // Explain type errors in more detail.
      "-feature", // Emit warning and location for usages of features that should be imported explicitly.
      "-unchecked", // Enable additional warnings where generated code depends on assumptions.
      "-Xcheckinit", // Wrap field accessors to throw an exception on uninitialized access.
      "-Xfatal-warnings", // Fail the compilation if there are any warnings.
      "-Xlint:adapted-args", // Warn if an argument list is modified to match the receiver.
      "-Xlint:constant", // Evaluation of a constant arithmetic expression results in an error.
      "-Xlint:delayedinit-select", // Selecting member of DelayedInit.
      "-Xlint:doc-detached", // A Scaladoc comment appears to be detached from its element.
      "-Xlint:inaccessible", // Warn about inaccessible types in method signatures.
      "-Xlint:infer-any", // Warn when a type argument is inferred to be `Any`.
      "-Xlint:missing-interpolator", // A string literal appears to be missing an interpolator id.
      "-Xlint:nullary-unit", // Warn when nullary methods return Unit.
      "-Xlint:option-implicit", // Option.apply used implicit view.
      "-Xlint:package-object-classes", // Class or object defined in package object.
      "-Xlint:poly-implicit-overload", // Parameterized overloaded implicit methods are not visible as view bounds.
      "-Xlint:private-shadow", // A private field (or class parameter) shadows a superclass field.
      "-Xlint:stars-align", // Pattern sequence wildcard must align with sequence component.
      "-Xlint:type-parameter-shadow", // A local type parameter shadows a type already in scope.
      "-Ywarn-dead-code", // Warn when dead code is identified.
      "-Ywarn-extra-implicit", // Warn when more than one implicit parameter section is defined.
      "-Ywarn-numeric-widen", // Warn when numerics are widened.
      "-Ywarn-unused:implicits", // Warn if an implicit parameter is unused.
      "-Ywarn-unused:imports", // Warn if an import selector is not referenced.
      "-Ywarn-unused:locals", // Warn if a local definition is unused.
      "-Ywarn-unused:params", // Warn if a value parameter is unused.
      "-Ywarn-unused:patvars", // Warn if a variable bound in a pattern is unused.
      "-Ywarn-unused:privates", // Warn if a private member is unused.
      "-Ybackend-parallelism", Math.min(Runtime.getRuntime.availableProcessors(), 8).toString, // Enable parallelization — scalac max is 16.
      "-Ycache-plugin-class-loader:last-modified", // Enables caching of classloaders for compiler plugins
      "-Ycache-macro-class-loader:last-modified", // and macro definitions. This can lead to performance improvements.
    )
  }
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
  def scalaVersion = "2.13.11"
  override def millSourcePath = super.millSourcePath / os.up
  override def mainClass = Some("com.fulcrumgenomics.sv.cmdline.SvMain")
  override def artifactName = "fgsv"
  def gitHash = Process("git rev-parse --short HEAD").lazyLines.head
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
    ivy"com.fulcrumgenomics:fgbio_2.13:2.2.1".excludeOrg(orgsToExclude:_*)
  )

  object test extends SbtModuleTests {
    override def ivyDeps = Agg(ivy"org.scalatest::scalatest:3.2.17")
    override def testFramework: Target[String] = T { "org.scalatest.tools.Framework" }

    // run mill tools.test.singleTest com.fulcrumgenomics.sv.x.y.x.TestClassName
    def singleTest(args: String*) = T.command {
      super.runMain("org.scalatest.run", args: _*)
    }
  }

  override def javacOptions = Seq("-source", "1.8", "-target", "1.8")

  def deployLocal = T { super.deployLocal(assembly(), "fgsv.jar")  }
}
