# If configured with --enable-mpi, compile mpi version of towhee
# and MPI compilers and libraries
# Otherwise, compile towhee and forcefield
bin_PROGRAMS = towhee forcefield
AM_FFLAGS = 

if UseMPI
MPI_CPPFLAGS = -DUSEMPI
MPI_LIBS = $(MPILIBS)
endif
# If configured with --enable-openkim, set the flag
if UseOpenkim
OK_CPPFLAGS = -DUSEOPENKIM
INCLUDES = -I${KIM_DIR}KIM_API
endif
# If configured with --enable-tramonto, set the flag
if UseTramonto
T_CPPFLAGS = -DUSETRAMONTO
endif
# If configure with --enable-lcao, set the flag
if UseLCAO
LCAO_CPPFLAGS = -DUSELCAO
endif
# If configure with --enable-internal, set the flag
if UseInternal
I_CPPFLAGS = -DUSEINTERNAL
endif
if UseSafeCompare
SC_CPPFLAGS = -DSAFE_COMPARE -Wall
SC_FFLAGS = -fbounds-check 
endif
if UseFixGNU
FG_FFLAGS = -fno-second-underscore
endif

# If configured with --disable-command-line-args, command line arguments
# won't be parsed by the towhee executable.
if DisableCLA
CLA_CPPFLAGS = -DDISABLE_CLA 
endif
#if UseDM
#DM_CPPFLAGS = -DUSE_DYNAMIC_MEMORY
#endif

AM_CPPFLAGS = $(MPI_CPPFLAGS) $(OK_CPPFLAGS) $(T_CPPFLAGS) $(LCAO_CPPFLAGS) $(I_CPPFLAGS)
AM_CPPFLAGS += $(SC_CPPFLAGS) $(CLA_CPPFLAGS) $(DM_CPPFLAGS)
AM_FFLAGS += $(SC_FFLAGS) $(FG_FFLAGS) $(DM_FFLAGS)
LIBS = $(MPI_LIBS)

# Add the Fortran libraries, so we can link with C
LDADD = @FLIBS@

towhee_SOURCES =  \
		accumulateaverages.F \
		allocate.F \
		angle_b_distribution.F \
		arbtocart.F \
		arccos.F \
		assemble.F \
		atd_cbmc_setup.F \
		atomshift.F \
		averages.F \
		bondorder.F \
		buildhelix.F \
		buildmolec.F \
		buildna.F \
		buildnanotube.F \
		buildpartial.F \
		buildprot.F \
		carttoarb.F \
		cbregrow.F \
		checkhmatrix.F \
		checkstruc.F \
		citations.F \
		compare.F \
		comparestruc.F \
		composite.F \
		cone.F \
		coneangle.F \
		conrot.F \
		createmolecule.F \
		crossproduct.F \
		ctrmas.F \
		dihedral_distribution.F \
		distance.F \
		dotproduct.F \
		energy.F \
		energybiasbox.F \
		energy_change.F \
		engatom.F \
		engtotal.F \
		expon.F \
		febias.F \
		fielddft.F \
		finalize.F \
		findrings.F \
		findtarget.F \
		gaussian.F \
		gaussprob.F \
		getcbangle.F \
		getcbbond.F \
		getcbdihed.F \
		getelement.F \
		getnbtype.F \
		getweight.F \
		globalaccess.F \
		globaldata.F \
		helixpos.F \
		initaverages.F \
		initconf.F \
		initialize.F \
		initloopclosure.F \
		integratedgaussprob.F \
		internal.F \
		interpolate.F \
		inverselawofcosine.F \
		inverthmatrix.F \
		kiss.F \
		lawofcosine.F \
		leastsquares.F \
		linclude.F \
		loopclosure.F \
		mainloop.F \
		mapconf.F \
		maxboxlength.F \
		mimage.F \
		minboxlength.F \
		monomers.F \
		openfile.F \
		pivot.F \
		polyamber96.F \
		polycharmm19.F \
		polycharmm22.F \
		polycharmm27.F \
		polycharmm27na.F \
		polyc19eef1.F \
		polyoplsaa.F \
		polysafetycheck.F \
		putarbinbox.F \
		putcartinbox.F \
		pupdate.F \
		qqassign.F \
		quantumenergy.F \
		radialpressure.F \
		ranlux.F \
		readclassical.F \
		readinput.F \
		readlammps.F \
		readquantum.F \
		readtowhee.F \
		readutils.F \
		recip.F \
		recippress.F \
		recipsum.F \
		resetcbmc.F \
		resetcoulomb.F \
		reseteamrho.F \
		rosenbluth.F \
		rotate.F \
		rotationmatrix.F \
		rwcbmc.F \
		rwconf.F \
		rwexternal.F \
		rwforcefield.F \
		rwsafe.F \
		scaled_potential.F \
		schedule.F \
		sdb_config.F \
		sdb_rdstruct.F \
		sdb_rdtoend.F \
		sdb_strget.F \
		sdb_strpars.F \
		setcbdihed.F \
		setclassical.F \
		setmixterms.F \
		setpotentiallogic.F \
		sphere.F \
		spline.F \
		stresstensor.F \
		swapmoves.F \
		switch.F \
		tail.F \
		tranatom.F \
		trancom.F \
		uniformbox.F \
		updatethreebondbias.F \
		updatetrmaxdisp.F \
		updatetwobondbias.F \
		updatevolmaxdisp.F \
		vangang.F \
		vangle.F \
		vbonbon.F \
		vbond.F \
		vcoulomb.F \
		veefone.F \
		vfield.F \
		vimproper.F \
		volnpt.F \
		volnvt.F \
		vonefive.F \
		vsasa.F \
		vembed.F \
		vthreebody.F \
		vtorsion.F \
		vtwobody.F \
		wcoulomb.F \
		writeangang.F \
		writeangle.F \
		writeaverages.F \
		writebond.F \
		writedlpoly.F \
		writedynamo.F \
		writeimproper.F \
		writeintro.F \
		writelammps.F \
		writelcao.F \
		writemovie.F \
		writenonbond.F \
		writeonefive.F \
		writepdb.F \
		writeruntime.F \
		writetorsion.F \
		writetowhee.F \
		writetramonto.F \
		wtwobody.F \
		control.c \
		dx-1597-2-7.c \
		globalinfo.c \
		jobfarm.c \
		mrg32k3a.c \
		rex.c \
		rnginit.c \
		scaled_potential_data.c \
		sturm.c \
		towhee.c

forcefield_SOURCES  = \
		arccos.F \
		averages.F \
		compare.F \
		crossproduct.F \
		distance.F \
		dotproduct.F \
		ffackl2004.F \
		ffalavi2005.F \
		ffamber96.F \
		ffaqvist.F \
		ffbelash2006.F \
		ffbelash2013.F \
		ffc19eef1.F \
		ffc19sasa.F \
		ffc27rigid.F \
		ffcatlowfaux.F \
		ffcharmm19.F \
		ffcharmm22.F \
		ffcharmm22fe.F \
		ffcharmm27.F \
		ffcharmm27x.F \
		ffcheck.F \
		ffclayff.F \
		ffcompass.F \
		ffcoon1987.F \
		ffcreate.F \
		ffcui1998.F \
		ffcui2002.F \
		ffdacnisua.F \
		ffdevel.F \
		ffdick1994.F \
		ffding1986.F \
		ffdreiding.F \
		ffdubb2004.F \
		ffdws2001.F \
		ffelli2002.F \
		ffepm.F \
		ffgala1994.F \
		ffgordon.F \
		ffgromos43a1.F \
		fffris2003.F \
		fffris2008.F \
		ffhardsphere.F \
                ffjaramillo.F \
		ffkfvbvs.F \
		fflast1993.F \
		fflgm.F \
		ffljium.F \
		ffmcy1976.F \
		ffmend2003.F \
		ffmgmstereo.F \
		ffmm2.F \
		ffmmff94.F \
		ffmorrow2002.F \
		ffnerdv1.F \
		ffnerdv2.F \
		ffnerdv3.F \
		ffoplsaa.F \
		ffoplsua.F \
		ffopls1996.F \
		ffopls2001.F \
		ffpana1989.F \
		ffpmf.F \
		ffpotter1997.F \
		ffqmff_viii.F \
		ffreadcharmmfile.F \
		ffreadsetflfile.F \
		ffrichar1995.F \
		ffshah2004.F \
		ffshukla1987.F \
		ffsks.F \
		ffsmith1994.F \
		ffsmmkmain.F \
		ffsmmknaip.F \
		ffspce.F \
		ffsquarewell.F \
		ffstil1985.F \
		ffsum2003.F \
		fftele1987.F \
		fftip3p.F \
		fftip4p.F \
		fftip5p.F \
		fftrappeeh.F \
		fftrappeua.F \
		fftrappeuaf.F \
		ffuff.F \
		ffunlu2004.F \
		ffvahid2010.F \
		ffvega1992.F \
		ffvink2001.F \
		ffwalt2001.F \
		ffweiner1984.F \
		ffweiner1986.F \
		ffwielop1985.F \
		getelement.F \
		getnbtype.F \
		globalaccess.F \
		globaldata.F \
		initialize.F \
		inverselawofcosine.F \
		inverthmatrix.F \
		openfile.F \
		readutils.F \
		rwforcefield.F \
		scaled_potential.F \
		globalinfo.c

EXTRA_towhee_SOURCES = \
		dvpm.h \
		functions.h \
		globalc.h \
		preproc.h \
		wrapperC.h \
		wrapperF.h 

###Fortran dependencies
accumulateaverages.o: preproc.h accumulateaverages.F
allocate.o: preproc.h allocate.F
angle_b_distribution.o: functions.h preproc.h angle_b_distribution.F
arbtocart.o: preproc.h arbtocart.F
arccos.o: functions.h arccos.F
assemble.o: preproc.h functions.h assemble.F
atd_cbmc_setup.o: preproc.h functions.h atd_cbmc_setup.F
atomshift.o: preproc.h functions.h atomshift.F
averages.o: preproc.h functions.h averages.F
bondorder.o: preproc.h functions.h bondorder.F
buildhelix.o: preproc.h functions.h buildhelix.F
buildmolec.o: preproc.h functions.h buildmolec.F
buildna.o: preproc.h functions.h buildna.F
buildnanotube.o: preproc.h functions.h buildnanotube.F
buildpartial.o: preproc.h buildpartial.F
buildprot.o: preproc.h functions.h buildprot.F
carttoarb.o: preproc.h carttoarb.F
cbregrow.o: preproc.h functions.h cbregrow.F
checkhmatrix.o: preproc.h functions.h checkhmatrix.F
checkstruc.o: preproc.h functions.h checkstruc.F 
citations.o: preproc.h citations.F
compare.o: preproc.h functions.h compare.F 
comparestruc.o: preproc.h functions.h comparestruc.F
composite.o: preproc.h composite.F
cone.o: preproc.h functions.h cone.F
coneangle.o: coneangle.F
conrot.o: preproc.h functions.h dvpm.h conrot.F
createmolecule.o: preproc.h functions.h createmolecule.F
crossproduct.o: preproc.h crossproduct.F
ctrmas.o: preproc.h functions.h ctrmas.F
dihedral_distribution.o: preproc.h functions.h dihedral_distribution.F
distance.o: functions.h distance.F
dotproduct.o: preproc.h functions.h dotproduct.F
energy.o: preproc.h functions.h energy.F
energybiasbox.o: preproc.h functions.h energybiasbox.F
energy_change.o: preproc.h functions.h energy_change.F
engatom.o: preproc.h functions.h engatom.F
engtotal.o: preproc.h functions.h engtotal.F
expon.o: preproc.h functions.h expon.F 
febias.o: preproc.h functions.h febias.F
ffackl2004.o: preproc.h functions.h ffackl2004.F
ffalavi2005.o: preproc.h functions.h ffalavi2005.F
ffamber96.o: preproc.h functions.h ffamber96.F
ffaqvist.o: preproc.h functions.h ffaqvist.F
ffbelash2006.o: preproc.h functions.h ffbelash2006.F
ffbelash2013.o: preproc.h functions.h ffbelash2013.F
ffc19eef1.o: preproc.h functions.h ffc19eef1.F
ffc19sasa.o: preproc.h functions.h ffc19sasa.F
ffc27rigid.o: preproc.h functions.h ffc27rigid.F
ffcatlowfaux.o: preproc.h functions.h ffcatlowfaux.F
ffcharmm19.o: preproc.h functions.h ffcharmm19.F
ffcharmm22.o: preproc.h functions.h ffcharmm22.F
ffcharmm22fe.o: preproc.h functions.h ffcharmm22fe.F
ffcharmm27.o: preproc.h functions.h ffcharmm27.F
ffcharmm27x.o: preproc.h ffcharmm27x.F
ffcheck.o: preproc.h ffcheck.F
ffclayff.o: preproc.h functions.h ffclayff.F
ffcompass.o: preproc.h functions.h ffcompass.F
ffcoon1987.o: preproc.h functions.h ffcoon1987.F
ffcreate.o: preproc.h functions.h ffcreate.F
ffcui1998.o: preproc.h functions.h ffcui1998.F
ffcui2002.o: preproc.h functions.h ffcui2002.F
ffdacnisua.o: preproc.h functions.h ffdacnisua.F
ffdevel.o: preproc.h functions.h ffdevel.F
ffdick1994.o: preproc.h functions.h ffdick1994.F
ffding1986.o: preproc.h functions.h ffding1986.F
ffdreiding.o: preproc.h functions.h ffdreiding.F
ffdubb2004.o: preproc.h functions.h ffdubb2004.F
ffdws2001.o: preproc.h functions.h ffdws2001.F
ffelli2002.o: preproc.h functions.h ffelli2002.F
ffepm.o: preproc.h functions.h ffepm.F
ffgala1994.o: preproc.h functions.h ffgala1994.F
ffgordon.o: preproc.h functions.h ffgordon.F
ffgromos43a1.o: preproc.h functions.h ffgromos43a1.F
fffris2003.o: preproc.h functions.h fffris2003.F
ffhardsphere.o: preproc.h functions.h ffhardsphere.F
ffjaramillo.o: preproc.h functions.h ffjaramillo.F
ffkfvbvs.o: preproc.h functions.h ffkfvbvs.F
fflast1993.o: preproc.h functions.h fflast1993.F
fflgm.o: preproc.h functions.h fflgm.F
ffljium.o: preproc.h functions.h ffljium.F
ffmcy1976.o: preproc.h functions.h ffmcy1976.F
ffmend2003.o: preproc.h functions.h ffmend2003.F
ffmgmstereo.o: preproc.h ffmgmstereo.F
ffmm2.o: preproc.h functions.h ffmm2.F
ffmmff94.o: preproc.h functions.h ffmmff94.F
ffmorrow2002.o: preproc.h functions.h ffmorrow2002.F
ffnerdv1.o: preproc.h functions.h ffnerdv1.F
ffnerdv2.o: preproc.h functions.h ffnerdv2.F
ffnerdv3.o: preproc.h functions.h ffnerdv3.F
ffoplsaa.o: preproc.h functions.h ffoplsaa.F
ffoplsua.o: preproc.h functions.h ffoplsua.F
ffopls1996.o: preproc.h functions.h ffopls1996.F
ffopls2001.o: preproc.h functions.h ffopls2001.F
ffpana1989.o: preproc.h functions.h ffpana1989.F
ffpmf.o: preproc.h functions.h ffpmf.F
ffpotter1997.o: preproc.h functions.h ffpotter1997.F
ffqmff_viii.o: preproc.h ffqmff_viii.F
ffreadcharmmfile.o: preproc.h functions.h ffreadcharmmfile.F
ffreadsetflfile.o: preproc.h functions.h ffreadsetflfile.F
ffshah2004.o: preproc.h functions.h ffshah2004.F
ffshukla1987.o: preproc.h functions.h ffshukla1987.F
ffsks.o: preproc.h functions.h ffsks.F
ffsmith1994.o: preproc.h functions.h ffsmith1994.F
ffsmmkmain.o: preproc.h functions.h ffsmmkmain.F
ffsmmknaip.o: preproc.h functions.h ffsmmknaip.F
ffspce.o: preproc.h functions.h ffspce.F
ffsquarewell.o: preproc.h functions.h ffsquarewell.F
ffstil1985.o: preproc.h functions.h ffstil1985.F
ffsum2003.o: preproc.h functions.h ffsum2003.F
fftele1987.o: preproc.h functions.h fftele1987.F
fftip3p.o: preproc.h functions.h fftip3p.F
fftip4p.o: preproc.h functions.h fftip4p.F
fftip5p.o: preproc.h functions.h fftip5p.F
fftrappeeh.o: preproc.h functions.h fftrappeeh.F
fftrappeua.o: preproc.h functions.h fftrappeua.F
fftrappeuaf.o: preproc.h functions.h fftrappeuaf.F
ffuff.o: preproc.h functions.h ffuff.F
ffunlu2004.o: preproc.h functions.h ffunlu2004.F
ffvahid2010.o: preproc.h functions.h ffvahid2010.F
ffvega1992.o: preproc.h functions.h ffvega1992.F
ffvink2001.o: preproc.h functions.h ffvink2001.F
ffwalt2001.o: preproc.h functions.h ffwalt2001.F
ffweiner1984.o: preproc.h functions.h ffweiner1984.F
ffweiner1986.o: preproc.h functions.h ffweiner1986.F
ffwielop1985.o: preproc.h functions.h ffwielop1985.F
fielddft.o: preproc.h fielddft.F
finalize.o: preproc.h finalize.F
findrings.o: preproc.h functions.h findrings.F
findtarget.o: preproc.h functions.h findtarget.F
gaussian.o: preproc.h functions.h gaussian.F
gaussprob.o: preproc.h functions.h gaussprob.F
getcbangle.o: preproc.h functions.h getcbangle.F
getcbbond.o: preproc.h functions.h getcbbond.F
getcbdihed.o: preproc.h functions.h getcbdihed.F
getelement.o: preproc.h getelement.F
getnbtype.o: preproc.h functions.h getnbtype.F
getweight.o: preproc.h functions.h getweight.F
globalaccess.o: preproc.h functions.h globalaccess.F
globaldata.o: preproc.h globaldata.F
helixpos.o: helixpos.F
initaverages.o: preproc.h initaverages.F
initconf.o: preproc.h initconf.F
initialize.o: preproc.h functions.h initialize.F
initloopclosure.o: preproc.h functions.h initloopclosure.F
integratedgaussprob.o: preproc.h functions.h integratedgaussprob.F
internal.o: preproc.h functions.h internal.F
interpolate.o: preproc.h functions.h interpolate.F
inverselawofcosine.o: functions.h inverselawofcosine.F
inverthmatrix.o: preproc.h inverthmatrix.F
lawofcosine.o: functions.h lawofcosine.F
leastsquares.o: leastsquares.F
linclude.o: preproc.h functions.h linclude.F 
loopclosure.o: preproc.h functions.h loopclosure.F
mainloop.o: preproc.h functions.h wrapperF.h mainloop.F
mapconf.o: preproc.h mapconf.F
maxboxlength.o: preproc.h functions.h maxboxlength.F
mimage.o: preproc.h functions.h mimage.F
minboxlength.o: preproc.h functions.h minboxlength.F
monomers.o: preproc.h monomers.F
openfile.o: preproc.h functions.h openfile.F
pivot.o: preproc.h dvpm.h functions.h pivot.F
polyamber96.o: preproc.h polyamber96.F
polycharmm19.o: preproc.h functions.h polycharmm19.F
polyc19eef1.o: preproc.h functions.h polyc19eef1.F
polycharmm22.o: preproc.h polycharmm22.F
polycharmm27.o: preproc.h polycharmm27.F
polycharmm27na.o: preproc.h polycharmm27na.F
polyoplsaa.o: preproc.h polyoplsaa.F
polysafetycheck.o: preproc.h polysafetycheck.F 
putarbinbox.o: preproc.h putarbinbox.F
putcartinbox.o: preproc.h putcartinbox.F
pupdate.o: preproc.h pupdate.F
qqassign.o: preproc.h functions.h qqassign.F
quantumenergy.o: preproc.h quantumenergy.F
radialpressure.o: preproc.h functions.h radialpressure.F
ranlux.o: ranlux.F
readclassical.o: preproc.h functions.h wrapperF.h readclassical.F
readinput.o: preproc.h wrapperF.h functions.h readinput.F
readlammps.o: preproc.h functions.h readlammps.F
readquantum.o: preproc.h functions.h readquantum.F
readtowhee.o: preproc.h wrapperF.h functions.h readtowhee.F
readutils.o: functions.h readutils.F
recip.o: preproc.h recip.F
recippress.o: preproc.h recippress.F
recipsum.o: preproc.h functions.h recipsum.F
resetcbmc.o: preproc.h resetcbmc.F
resetcoulomb.o: preproc.h functions.h resetcoulomb.F
reseteamrho.o: preproc.h functions.h reseteamrho.F
rex.o: globalc.h preproc.h rex.c
rosenbluth.o: preproc.h functions.h rosenbluth.F
rotate.o: preproc.h functions.h rotate.F
rotationmatrix.o: preproc.h rotationmatrix.F
rwcbmc.o: preproc.h functions.h rwcbmc.F
rwconf.o: preproc.h functions.h rwconf.F
rwexternal.o: preproc.h rwexternal.F
rwforcefield.o: preproc.h functions.h rwforcefield.F
rwsafe.o: preproc.h functions.h rwsafe.F
scaled_potential.o: preproc.h scaled_potential.F
schedule.o: preproc.h functions.h schedule.F
sdb_config.o: sdb_config.F
sdb_rdstruct.o: sdb_rdstruct.F
sdb_rdtoend.o: sdb_rdtoend.F 
sdb_strget.o: sdb_strget.F
sdb_strpars.o: sdb_strpars.F
setcbdihed.o: preproc.h functions.h setcbdihed.F
setclassical.o: preproc.h functions.h setclassical.F
setmixterms.o: preproc.h functions.h setmixterms.F
setpotentiallogic.o: preproc.h setpotentiallogic.F
sphere.o: functions.h sphere.F
spline.o: preproc.h functions.h spline.F
stresstensor.o: preproc.h functions.h stresstensor.F
swapmoves.o: preproc.h functions.h swapmoves.F
switch.o: preproc.h switch.F
tail.o: preproc.h functions.h tail.F
tranatom.o: preproc.h dvpm.h functions.h tranatom.F
trancom.o: preproc.h functions.h trancom.F
uniformbox.o: functions.h uniformbox.F
updatemaxtrdisp.o: preproc.h updatemaxtrdisp.F
updatethreebondbias.o: preproc.h functions.h updatethreebondbias.F
updatetwobondbias.o: preproc.h functions.h updatetwobondbias.F
updatevolmaxdisp.o: functions.h preproc.h updatevolmaxdisp.F
vangang.o: preproc.h functions.h vangang.F
vangle.o: preproc.h functions.h vangle.F
vbonbon.o: preproc.h functions.h vbonbon.F
vbond.o: preproc.h functions.h vbond.F
vcoulomb.o: preproc.h functions.h vcoulomb.F
veefone.o: preproc.h functions.h veefone.F
vfield.o: preproc.h functions.h vfield.F
vimproper.o: preproc.h functions.h vimproper.F
volnpt.o: preproc.h functions.h volnpt.F
volnvt.o: preproc.h functions.h volnvt.F
vonefive.o: preproc.h functions.h vonefive.F
vsasa.o: preproc.h functions.h vsasa.F
vembed.o: preproc.h functions.h vembed.F
vthreebody.o: preproc.h functions.h vthreebody.F
vtorsion.o: preproc.h functions.h vtorsion.F
vtwobody.o: preproc.h functions.h vtwobody.F
wcoulomb.o: preproc.h functions.h wtwobody.F
writeangang.o: preproc.h functions.h writeangang.F
writeangle.o: preproc.h functions.h writeangle.F
writeaverages.o: preproc.h writeaverages.F
writebond.o: preproc.h functions.h writebond.F
writedlpoly.o: preproc.h functions.h writedlpoly.F
writedynamo.o: preproc.h functions.h writedynamo.F
writeimproper.o: preproc.h functions.h writeimproper.F 
writeintro.o: preproc.h writeintro.F 
writelammps.o: preproc.h functions.h writelammps.F
writelcao.o: preproc.h functions.h writelcao.F
writemovie.o: preproc.h functions.h writemovie.F
writenonbond.o: preproc.h functions.h writenonbond.F
writeonefive.o: preproc.h functions.h writeonefive.F
writepdb.o: preproc.h functions.h writepdb.F
writeruntime.o: preproc.h writeruntime.F
writetorsion.o: preproc.h functions.h writetorsion.F
writetowhee.o: preproc.h functions.h writetowhee.F
writetramonto.o: preproc.h functions.h writetramonto.F
wtwobody.o: preproc.h functions.h wtwobody.F
###C dependencies
control.o: globalc.h control.c
dx-1597-2-7.o: globalc.h dx-1597-2-7.c
globalinfo.o: globalc.h preproc.h globalinfo.c
jobfarm.o: globalc.h jobfarm.c
mrg32k3a.o: globalc.h mrg32k3a.c
rex.o: globalc.h rex.c
rnginit.o: globalc.h preproc.h rnginit.c
scaled_potential_data.o: globalc.h scaled_potential_data.c
sturm.o: sturm.c
towhee.o: globalc.h towhee.c wrapperC.h
