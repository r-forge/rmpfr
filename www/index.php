
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="/"><img src="<?php echo $themeroot; ?>/images/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<h2>Rmpfr: Multiple Precision Floating-point in R</h2>
<p>Multiple precision numbers &amp; computations for
 <a href="http://www.r-project.org">R</a> via GNU MPFR
 <a href="http://www.mpfr.org">www.mpfr.org</a>
 and <a href="http://gmplib.org">GMP</a>.
</p>

<p>
<samp>Rmpfr</samp> provides S4 classes and methods for
arithmetic including transcendental ("special") functions for
arbitrary precision floating point numbers.
To this end, it interfaces to
the LGPL ed <a href="http://www.mpfr.org">MPFR </a>(Multiple Precision Floating-Point Reliable) Library which itself is based on the  <a href="http://gmplib.org">GMP</a> (GNU Multiple Precision) Library.
</p>

<p> The (Swiss mirror of the)
 <a href="http://cran.r-project.org/">CRAN</a> page of <samp>Rmpfr</samp>,
 <a href="http://stat.ethz.ch/CRAN/web/packages/Rmpfr/">Rmpfr CRAN page</a>
    contains the reference manual (help pages of all functions), vignettes
    (in pdf) and further relevant information about the package.
</p>

<p>
  In 2012, Rmpfr is also interfacing to the
  <a href="http://mulcyber.toulouse.inra.fr/projects/gmp">
  R package <samp>gmp</samp></a>.  gmp implements <strong>exact</strong>
  arithmetic for large integers and rationals, the former also <em>"modulo m"</em>,
<!--  i.e., (mathematically) <em>"in Z / mZ"</em>. -->
</p>

<!-- end of project description -->

<p> The <samp>Rmpfr</samp> <strong>project summary page</strong> is <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>

</body>
</html>
