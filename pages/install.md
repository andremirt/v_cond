---
layout: page
title: Install
permalink: install/
---

To compute \\(v_{cond}(\mathbf{r})\\) you will have to install

<ul>
  <li>our modified vesion of GAMESS-US</li>
  <li>the SLATEC library</li>
  <li>our v_cond routines</li>
</ul>

Find below a Linux sample installation for the <a href="https://surfsara.nl/systems/lisa">Lisa cluster</a> of the dutch high performance computing center SURFsara.

<h2>Prerequisites</h2>

You will need

<ul>
  <li>a Fortran compiler</li>
  <li>a C compiler</li>
  <li>the <a href="http://www.netlib.org/blas/">BLAS routines</a>, e.g. from <a  href="http://math-atlas.sourceforge.net/">ATLAS</a> or <a href="https://software.intel.com/en-us/intel-mkl">Intel's MKL</a></li>
</ul>

For our Lisa installation make them available with the command

<pre>
<code class="lb">module load fortran/intel c/intel/64 mkl/64</code>
</pre>

<h2>GAMESS-US Install</h2>

Please register as GAMESS-US user online at the <a href="http://www.msg.ameslab.gov/gamess/License_Agreement.html">Gordon research group at Iowa State University, US</a>. You can then obtain the most recent version of GAMESS-US. For our purposes we will install a modiefied version of GAMESS-US that builds up on the GAMESS-US version 1 OCT 2010 (R1).

Extensive installation instructions for GAMESS-US can be found in the <a href="https://github.com/andremirt/v_cond/blob/master/gamess/misc/readme.unix">./gamess/misc/readme.unix</a>.

Here we continue with the Lisa installation. 

<h5>Configure GAMESS-US installation</h5>

To configure your installation run
<pre>
<code  class="language-bash" data-lang="bash">cd gamess
csh config</code>
</pre>

For Lisa we set
<ul>
<table>
  <tr>
    <td>target machine:</td>
    <td>linux64</td>
  </tr>
  <tr>
    <td>Fortran:</td>
    <td>ifort</td>
  </tr>
  <tr>
    <td>version:</td>
    <td>13</td>
  </tr>
  <tr>
    <td>math libraries:</td>
    <td>mkl</td>
  </tr>
  <tr>
    <td>mkl path:</td>
    <td>/sara/sw/mkl-11.0.2/composer_xe_2013.2.146/mkl/</td>
  </tr>
  <tr>
    <td>communication library:</td>
    <td>sockets</td>
  </tr>
</table>
</ul>

<h5>Compile GAMESS-US distributed data interface (ddi)</h5>

To compile the ddi run
<pre>
<code  class="language-bash" data-lang="bash">cd ddi
csh compddi</code>
</pre>
Then move <code>ddikick.x</code> to the GAMESS-US installation directory
<pre>
<code  class="language-bash" data-lang="bash">mv ddikick.x ..</code>
</pre>
and leave the ddi directory
<pre>
<code  class="language-bash" data-lang="bash">cd ..</code>
</pre>

<h5>Compile and link GAMESS-US</h5>

To compile the GAMESS-US routines run
<pre>
<code  class="language-bash" data-lang="bash">csh compall</code>
</pre>

If needed, add execution permision to <code>comp</code>
<pre>
<code  class="language-bash" data-lang="bash">chmod u+x comp</code>
</pre>

Before linking GAMESS-US the <code>install.info</code> needs to be modified
change 
<pre>
<code  class="language-bash" data-lang="bash">setenv GMS_MATHLIB_PATH    /sara/sw/mkl-11.0.2/composer_xe_2013.2.146/mkl///lib/em64t</code>
</pre>
to
<pre>
<code  class="language-bash" data-lang="bash">setenv GMS_MATHLIB_PATH    /sara/sw/mkl-11.0.2/composer_xe_2013.2.146/mkl/lib/intel64</code>
</pre>
and add the line
<pre>
<code  class="language-bash" data-lang="bash">setenv GMS_MKL_VERNO       11</code>
</pre>
Finally, link GAMESS-US
<pre>
<code  class="language-bash" data-lang="bash">csh lked gamess</code>
</pre>

<h5>Modify execution script</h5>

In the file <code>rungmsvcond</code> set the <code>SCR</code> and <code>PATH_TO_GAMESS</code> variables.

A choice for the Lisa <code>SCR</code> is
<pre>
<code  class="language-bash" data-lang="bash">set SCR=/home/$USER/scr/</code>
</pre>

where you need to make sure that the directory <code>/home/$USER/scr/</code> exists.

Add execution permision to <code>rungmsvcond</code>
<pre>
<code  class="language-bash" data-lang="bash">chmod u+x rungmsvcond</code>
</pre>

And, if desired, create a link to the binary file in <code>~/bin</code>
<pre>
<code  class="language-bash" data-lang="bash">ln -s rungmsvcond ~/bin/.</code>
</pre>

<h2>Compiling the SLATEC library</h2>
<a href="www.netlib.org/slatec/">SLATEC Common Mathematical Library</a> is a comprehensive software library containing over 1400 general purpose mathematical and statistical routines written in Fortran 77. Here we compile the static library of a slightly modiefied version that is adapted to the Lisa cluster.

Enter the SLATEC static directory
<pre>
<code  class="language-bash" data-lang="bash">cd ../slatec/static</code>
</pre>

Type 
<pre>
<code  class="language-bash" data-lang="bash">make</code>
</pre>
to compile the routines and create the library.

<h2>Compiling the v_cond routines</h2>

Enter the v_cond directory
<pre>
<code  class="language-bash" data-lang="bash">cd ../../v_cond</code>
</pre>
and run
<pre>
<code  class="language-bash" data-lang="bash">. comp.sh</code>
</pre>
to compile <code>vcond.exe</code>. Eventually create a link to the binary file in <code>~/bin</code>
<pre>
<code  class="language-bash" data-lang="bash">ln -s vcond.exe ~/bin/.</code>
</pre>

