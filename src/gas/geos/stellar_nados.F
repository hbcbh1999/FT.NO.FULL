#include <gdecs/fortdecs.h>
      subroutine nadstate(ftden,fttemp,ftpres,ftener,ftentr,ftxmass,
     1 ftaion, ftzion,ftion)      
      implicit none
      save
      include 'vector_eos.dek'

c..
c. .tests the nadyozhin eos routine
c..
c..ionmax  = number of isotopes in the network
c..xmass   = mass fractions
c..ymass   = molar fractions
c..aion    = number of nucleons
c..zion    = number of protons

      integer          loopval,ionmaxsize
      parameter        (ionmaxsize=20)
      double precision xmass(ionmaxsize),ymass(ionmaxsize),
      1                 aion(ionmaxsize),zion(ionmaxsize),temp,den,
      2                 abar,zbar

      double precision      ftden,fttemp,ftpres,ftener,ftentr,
      1                       ftxmass(*),ftaion(*), ftzion(*)
      integer ftion


cc      real                  ftxmass(*),ftaion(*), ftzion(*)
c      write(*,*)  'inside nadstate'

c..set the mass fractions, z's and a's of the composition
c..hydrogen
c      xmass(1) = ftxmass(1) 
c      aion(1)  =  
c      zion(1)  = 1.0d0

c..helium
c      xmass(2) = ftxmass(2) 
c      aion(2)  = 4.0d0
c      zion(2)  = 2.0d0

c..carbon 12
c      xmass(3) = ftxmass(3) 
c      aion(3)  = 12.0d0
c      zion(3)  = 6.0d0
c       write(*,*) ftion
       do loopval = 0, ftion, 1
       xmass(loopval) = ftxmass(loopval)
       aion(loopval)  = ftaion(loopval) 
       zion(loopval)  = ftzion(loopval)
       end do 

c      write(*,*)  'inside nadstate'
c..get abar, zbar and a few other composition variables
      call azbar(xmass,aion,zion,ftion,
      1           ymass,abar,zbar)


c      write(*,*)  'inside nadstate'
c..set the input vector
      temp_row(1) = fttemp 
      den_row(1)  = ftden 
      abar_row(1) = abar
      zbar_row(1) = zbar
c..here the pipeline is only 1 element long
      jlo_eos = 1
      jhi_eos = 1


c..call the eos

c      write(*,*)  'inside nadstate'
      call nados

      ftpres = ptot_row(1)
      ftener = etot_row(1)
      ftentr = stot_row(1)
c..write out the results
c      call pretty_eos_out('nados: ')

c      stop 'normal termination'
      end   






      subroutine azbar(xmass,aion,zion,ionmax,
     1                 ymass,abar,zbar)
      implicit none
      save

c..this routine calculates composition variables for an eos routine

c..input:
c..mass fractions     = xmass(1:ionmax)
c..number of nucleons = aion(1:ionmax)
c..charge of nucleus  = zion(1:ionmax)
c..number of isotopes = ionmax

c..output:
c..molar abundances        = ymass(1:ionmax), 
c..mean number of nucleons = abar
c..mean nucleon charge     = zbar


c..declare
      integer          i,ionmax,ionmaxsize
      parameter        (ionmaxsize=20)
      double precision xmass(ionmaxsize),aion(ionmaxsize),
     1                 ymass(ionmaxsize),abar,zbar,zbarxx,ytot1,
     2                 zion(ionmaxsize)

      zbarxx  = 0.0d0
      ytot1   = 0.0d0
      do i=1,ionmax
       ymass(i) = xmass(i)/aion(i)
       ytot1    = ytot1 + ymass(i)
       zbarxx   = zbarxx + zion(i) * ymass(i)
      enddo
      abar   = 1.0d0/ytot1
      zbar   = zbarxx * abar
      return
      end






      subroutine pretty_eos_out(whose)
      implicit none
      save
      include 'vector_eos.dek'
c..
c..writes a pretty output for the eos tester
c..
c..declare
      integer     j
      character*7 whose

c..popular formats
01    format(1x,t2,a,t11,'total',t24,'ion',t34,'e- & e+',
     1       t46,'radiation',t58,'coulomb')
02    format(1x,t2,a,1p6e12.4)
03    format(1x,t2,a6,1pe12.4,t22,a6,1pe12.4,
     1         t42,a6,1pe12.4,t62,a6,1pe12.4)



      do j=jlo_eos,jhi_eos


c..the input 
      write(6,03) 'temp =',temp_row(1),'den  =',den_row(1),
     1            'abar =',abar_row(1),'zbar =',zbar_row(1)
      write(6,*) ' ' 


c..and the output
c..first the totals from each of the components
      write(6,01)  whose
      write(6,02) 'pres =',
     1            ptot_row(j),pion_row(j),pele_row(j),
     2            prad_row(j),pcou_row(j)
      write(6,02) 'ener =',
     1            etot_row(j),eion_row(j),eele_row(j),
     2            erad_row(j),ecou_row(j)
      write(6,02) 'entr =',
     1            stot_row(j),sion_row(j),sele_row(j),
     2            srad_row(j),scou_row(j)

c..derivatives of the totals with respect to the input variables
      write(6,*)  ' '
      write(6,03) 'dp/dd=',dpd_row(j),'dp/dt=',dpt_row(j),
     1            'dp/da=',dpa_row(j),'dp/dz=',dpz_row(j)
      write(6,03) 'de/dd=',ded_row(j),'de/dt=',det_row(j),
     1            'de/da=',dea_row(j),'de/dz=',dez_row(j)
      write(6,03) 'ds/dd=',dsd_row(j),'ds/dt=',dst_row(j),
     1            'ds/da=',dsa_row(j),'ds/dz=',dsz_row(j)


c..derivatives of the electron-positron compoenets with
c..respect to the input variables
      write(6,*) ' ' 
      write(6,03) 'dpepd=',dpepd_row(j),'dpept=',dpept_row(j),
     1            'dpepa=',dpepa_row(j),'dpepz=',dpepz_row(j)
      write(6,03) 'deepd=',deepd_row(j),'deept=',deept_row(j),
     1            'deepa=',deepa_row(j),'deepz=',deepz_row(j)
      write(6,03) 'dsepd=',dsepd_row(j),'dsept=',dsept_row(j),
     1            'dsepa=',dsepa_row(j),'dsepz=',dsepz_row(j)


c..the thermodynamic consistency relations, these should all be
c..at the floating poiint limit of zero
      write(6,*) ' ' 
      write(6,03) 'maxw1=',dse_row(j),'maxw2=',dpe_row(j),
     1            'maxw3=',dsp_row(j)


c..number density of electrons, poistrons, matter electrons, and ions
      write(6,03) 'xne  =',xne_row(j),'xnp  =',xnp_row(j),
     1            'xnem =',xnem_row(j),'xni  =',xni_row(j)


c..derivatibves of the electron number density with 
c..respect to the input variables
      write(6,03) 'dxned=',dxned_row(j),'dxnet=',dxnet_row(j),
     1            'dxnea=',dxnea_row(j),'dxnez=',dxnez_row(j)


c..electron chemical potential, positron chemical potential
c..and derivatives of electron chemical potential with respect
c..to the input variables
      write(6,03) 'eta  =',etaele_row(j),'etap =',etapos_row(j)
      write(6,03) 'detad=',detad_row(j),'detat=',detat_row(j),
     1            'detaa=',detaa_row(j),'detaz=',detaz_row(j)


c..specific heats, and ratio of electostatic to thermal energy
      write(6,03) 'cp   =',cp_row(j),'cv   =',cv_row(j),
     1            'plasg=',plasg_row(j)

c..the 3 gammas and the sound speed
      write(6,03) 'gam1 =',gam1_row(j),'gam2 =',gam2_row(j),
     1            'gam3 =',gam3_row(j),'csond=',cs_row(j)
      write(6,*) ' '

      enddo
      return
      end






      subroutine nados

c*********************************************************************** 
c  *** equation of state for completely ionized matter 
c  *** electron & positron component --- fermi-dirac statistics using 
c                       various asymptotic expansions where possible 
c  *** ion component --- a perfect gas approximation 
c  *** black-body radiation 
c*********************************************************************** 
c                            references 
c   1. nadyozhin d.k. 1974, "naucnye informatsii", nos. 32, 33 
c                            (ucsc science library: qb 1 a4) 
c*********************************************************************** 

      implicit double precision (a-h,o-z) 
      save
      include 'vector_eos.dek'


      double precision chit,chid,cv,cp,gam1,gam2,gam3,nabad,sound,zz,z,
     1                 prad,pion,dpiondd,dpiondt,erad,eion,srad,sion
      double precision clight,con2,me,mecc,kerg,avo
      parameter        (clight = 2.997925d10,
     1                  con2 = clight * clight,
     2                  me   = 9.1093897d-28, 
     3                  mecc = me * con2,
     4                  kerg    = 1.380658d-16,
     5                  avo     = 6.0221367d23)


c 
c  *** the arguments 
      common/arg/t,den,psi 
      equivalence(den,pl) 
      common/iarg/lst,kentr,kpar,jurs,jkk 
c*********************************************************************** 
c 
c  ***   t --- temperature in 10**9 k 
c  *** den --- density in 10**7 g/ccm 
c  *** psi --- parameter of degeneracy. works as an argument only when 
c              one enters entry fd12f to get fermi-dirac functions, 
c              otherwise it is calculated as a function of t and den. 
c  *** lst, kentr, kpar --- the keys aimed to make the calculations 
c               faster where possible (default: lst, kentr, kpar = 0) 
c  *** lst=0 --- epeos calculates only thermodynamic functions p,e,s 
c      lst=1 --- the first temperature-derivatives are calculated 
c                in addition to the thermodynamic functions 
c      lst=2 --- extends calculations to get density-derivatives as well 
c  *** kentr=0 --- turns the calculation of entropy off and suspends 
c            the calculation of psi in a perfect gas asymptotics (nz=1). 
c      kentr=1 --- turns the calculation of entropy on. 
c  *** kpar=0 --- when in relativistic asymptotics (nz=4), turns off the 
c        calculation of total number of pair (hpr),(kpar=1 --- turns on) 
c      kpar is inactive for other asymptotics. 
c 
c  *** jkk --- the current number of mesh point, is inactive in this  
c        version of epeos. however, it appears in epeos error messages 
c        and, if defined, may be useful for locating of errors in 
c        a program calling epeos. 
c*********************************************************************** 
      common/nz/nz 
c*********************************************************************** 
c  *** nz --- specifies the operational mode (default: nz=0) 
c      nz=0 --- calculations with the overall search for 
c               the appropriate working region 
c      for 0<nz<6 epeos works within one of five following modes 
c      independent of the values of temperature and density specified 
c      nz=1 --- perfect gas approximation with the first order 
c               corrections for degeneracy and pairs 
c      nz=2 --- expansion over half-integer fermi--dirac functions 
c      nz=3 --- chandrasekhar's expansion for degenerate gas 
c      nz=4 --- relativistic asymptotics 
c      nz=5 --- quadratures taken with the gauss method 
c*********************************************************************** 
      common/az/as,zs,scn 
c*********************************************************************** 
c  ***  as --- mean mass number, zs --- mean nuclear charge 
c          emue=as/zs --- mean electron molecular weight: the total 
c          number of "atomic" electrons in a unit volume, nea, 
c          amounts to (density)/(mu*emue), mu is atomic mass unit. 
c          for a mixture:   as=1/sum{y_i},  zs=as*sum{z_i*y_i}, where 
c          y_i=x_i/a_i,  and  x_i  being a mass fraction of 'i' species. 
c  *** scn --- additive entropy constant for the ion component 
c          for a gas of nuclei with mass number a, and spin i 
c          scn=ln[(2i+1)*a**2.5], for iron-56, scn=2.5*ln(56)=10.063379. 
c          for a mixture:   scn=as*sum{y_i*ln[(2i_i+1)*(a_i)**2.5}. 
c  *** jurs --- if jurs=0 then common-block  /az/ is to be preliminary 
c            filled with all necessary information. (default values 
c            of as,zs,scn are specified in block data for pure iron-56). 
c if jurs=1, epeos applies to subroutine chemic for as,zs,scn 
c*********************************************************************** 
c 
c                         diagnostics 
c*********************************************************************** 
c epeos opens file 'epeos.mes', writes the epeos version in it, analyses 
c the values of arguments in common-blocks /arg/, /iarg/, /nz/, /az/ and 
c then writes additional information in 'epeos.mes' if something wrong 
c with the arguments: in particular, epeos stops when  t = or < 0. 
c in case  den < 0 , epeos sends a warning in 'epeos.mes' and continues 
c to calculate with den=0 (black body radiation only). 
c*********************************************************************** 
c 
c  *** the results of calculations 
      common/result/p,e,s,sk,pt,et,st 
      common/str/ppl,epl,spl,cp,gam,da,dpe,dse,dsp,beta 
      common/resel/pe,ee,se,sek,hpr 
      common/nzr/nzr 
c*********************************************************************** 
c  *** p --- total pressure in units 10**24 erg/ccm 
c  *** e --- total specific energy in units 10**17 erg/g 
c  *** s --- total entropy in units  10**8 erg/(g k) 
c  *** sk --- total dimensionless entropy per nucleon: sk=s*mu/kb 
c             mu -- atomic mass unit, kb -- boltzmann's constant 
c  *** pt,et,st --- temperature derivatives of p,e,s at constant density 
c  *** ppl,spl --- density derivatives of p and s at constant temperature 
c  *** epl---- density derivatives of e multiplied by density den 
c  *** pe,ee,se,sek ----  pressure, specific energy and entropy 
c                         of the electron-positron gas component 
c  *** hpr --- the total number of the electron-positron pairs 
c              per "atomic" electron (hpr=mu*emue*np/den), where 
c              'np' is the number of pairs per unit volume.  
c  *** cp --- specific heat at constant pressure 
c  *** gam --- adiabatic index = {d log(p)/d log(den)} at s=const 
c  *** da --- logarithmic adiabatic temperature gradient 
c  *** dpe = (den/p)(epl+t*pt/den)-1 -- thermodynamic identity: dpe=0 
c  *** dse = t*st/et-1 -- thermodynamic identity: dse=0 
c  *** dsp = -spl*(den/pt)-1 -- thermodynamic identity: dsp=0 
c  *** beta --- ratio of black-body radiation pressure to total pressure 
c  *** nzr --- identificator of working region on t-den plane when nz=0 
c*********************************************************************** 
      common/fdf/f12,f32,f52,f72,f12s,f32s,f52s,f72s 
c*********************************************************************** 
c  *** f12,f32,f52,f72 --- half-integer fermi-dirac functions 
c  *** f12s,f32s,f52s,f72s --- the first derivatives of f-d functions 
c  *** psi (in common/arg/t,den,psi) is the argument of f-d functions 
c      there exists special entry to get these functions separately -- 
c      use command call fd12f after specifying psi in common/arg/ 
c*********************************************************************** 
      dimension cu(62),ck(3),a(17),c(8),b(22),b1(4),b2(4),b3(4),b4(4), 
     1b5(4),c1(4),c2(4),c3(4),c4(4),c5(4),c6(4),d1(4),d2(4),d3(4),d4(4), 
     2d5(4),d6(4),d(4),a1(4),a2(4),a3(4),a4(4),df1(4),df2(4) 
      dimension uio(5),ui1(5),ui2(5),cio(5),ci1(5),ci2(5),aio(5),ai1(5) 
      dimension ai2(5),xxi(5),aai(5),cci(5),bbi(5),wk1(5),wk2(5),wk3(5) 
      dimension uac(95),wk4(5),wk5(5),wk6(5) 
      dimension cpp(5),abc(85),ado(5),ad1(5),ad2(5),bdo(5),bd1(5),fgs(8) 
      dimension bd2(5),cdo(5),cd1(5),cd2(5),gdo(5),gd1(5),gd2(5) 
      dimension ggsi(5),zzi(5),vvi(5),hhi(5),ggi(5) 
      dimension asp(3),bsp(3),csp(3),gsp(3),aspa(3),bspa(3),cspa(3),gspa 
     1(3),abcg(24),wk7(5) 
      equivalence(psi,pc) 
      equivalence (uio(1),uac(1)),(ui1(1),uac(6)),(ui2(1),uac(11)), 
     1(cio(1),uac(16)),(ci1(1),uac(21)),(ci2(1),uac(26)),(aio(1),uac(31) 
     2),(ai1(1),uac(36)),(ai2(1),uac(41)),(xxi(1),uac(46)),(aai(1),uac(5 
     31)),(cci(1),uac(56)),(bbi(1),uac(61)),(wk1(1),uac(66)),(wk2(1),uac 
     4(71)),(wk3(1),uac(76)),(wk4(1),uac(81)),(wk5(1),uac(86)),(wk6(1),u 
     5ac(91)) 
      equivalence(abc(1),ado(1)),(abc(6),ad1(1)),(abc(11),ad2(1)), 
     1(abc(16),bdo(1)),(abc(21),bd1(1)),(abc(26),bd2(1)),(abc(31),cdo(5) 
     2),(abc(36),cd1(1)),(abc(41),cd2(1)),(abc(46),gdo(1)),(abc(51),gd1( 
     31)),(abc(56),gd2(1)),(abc(61),ggsi(1)),(abc(66),zzi(1)), 
     4(abc(71),vvi(1)),(abc(76),hhi(1)),(abc(81),ggi(1)),(abcg(1),asp(1) 
     5),(abcg(4),bsp(1)),(abcg(7),csp(1)),(abcg(10),gsp(1)),(abcg(13),as 
     6pa(1)),(abcg(16),bspa(1)),(abcg(19),cspa(1)),(abcg(22),gspa(1)) 


c..for the coulomb corrections
      double precision dsdd,dsda,lami,inv_lami,lamidd,ktinv,
     1                 plasg,plasgdd,plasgdt,
     2                 aa1,bb1,cc1,dd1,ee1,aa2,bb2,cc2,
     3                 ecoul,decouldd,decouldt,
     4                 pcoul,dpcouldd,dpcouldt,
     5                 scoul,dscouldd,dscouldt,
     6                 third,qe,esqu,forth,pie
      parameter        (third = 1.0d0/3.0d0,
     1                  qe    = 4.8032068d-10,
     2                  esqu  = qe * qe,
     3                  forth = 4.0d0/3.0d0,
     4                  pie   = 3.1415926535897932384d0)

c..see yakovlev & shalybkov 1989
      data aa1,bb1,cc1,dd1,ee1,aa2,bb2,cc2
     1     /-0.898004d0, 0.96786d0, 0.220703d0, -0.86097d0,
     2       2.5269d0  , 0.29561d0, 1.9885d0,    0.288675d0/



      data eit/1.d-6/,dst/1.d-3/,grcf/1.3d0/,grif/0.7d0/,tpar/0.3d0/, 
     1tt1/0.8d0/,tt2/0.96d0/,tpar1/0.32d0/,ro1/0.12d0/,ro2/3.4884680d-2/, 
     2ep2/0.7d0/,ps1/7.2427389d-5/,ro1l/-2.1202635d0/,ro2l/-3.3557075d0/, 
     3xep2/0.3d0/,rotp/3.1220424d-10/,pss1/0.1662467d0/, 
     4 pss2/0.1881164d0/ 
      data cu/5.93013d0,0.831434d0,1.5d0,1.d0,2.8125d0,5.5957031d0, 
     13.75d0,5.15625d0,1.25d0,0.9375d0,1.8652344d0,5.25d0,1.265625d1, 
     20.5625d0,2.d0,0.5d0,3.d0,2.25d0,0.75d0,2.984375d0,1.3125d1, 
     31.6875d1,3.5d0,2.5d0,1.6787109d1,1.40625d0,2.5d0,0.25d0, 
     45.6568542d0,5.75d0,1.453125d1,4.63125d1,5.8125d1,1.35d1, 
     51.013429d-3,1.33333333333d0,4.d0,2.8613184d-2,0.94051952d-1, 
     60.71428571d0,0.39650038026d-1,0.21875d0,0.66666666667d0,0.6d0, 
     71.05d0,0.18007375d0,0.292178d0,0.3601475d0,0.33333333333d0, 
     85.d0,8.d0,2.66666666667d0,24.d0,15.d0,6.d0,9.d0,1.2d1, 
     91.02677135171d+1,4.9305117064d0,4.80195683116d-1, 
     a8.09755744169d-2,7.99196885645d-1/ 
      data a1/6.78196814d-1,6.78183752d-1,6.94779549d-1,7.60042563d-1/ 
      data a2/5.37738527d-1,5.37346659d-1,4.87559266d-1,4.09243650d-1/ 
      data a3/1.73981666d-1,1.70062988d-1,2.19850381d-1,0.251176627d0/ 
      data a4/2.38231503d-2,1.07608902d-2,-5.83490747d-3,-.0100117403d0/ 
      data df1/3.47963332d-1,3.40125976d-1,4.39700762d-1,0.502353254d0/ 
      data df2/7.14694509d-2,3.22826706d-2,-1.75047241d-2, 
     d         -3.00352209d-2/ 
      data  b1/1.15292705d0,1.15292656d0,1.14670313d0,1.08551906d0/ 
      data  b2/1.01729522d0,1.01727563d0,1.04216932d0,1.14006384d0/ 
      data  b3/4.03303895d-1,4.03009994d-1,3.65669449d-1,3.06932737d-1/ 
      data b4/8.69908330d-2,8.50314940d-2,1.09925190d-1,1.25588313d-1/ 
      data b5/8.93368136d-3,4.03533382d-3,-2.18809030d-3,-3.75440261d-3/ 
      data c1/3.08286343d0,3.08286341d0,3.08597512d0,3.16245521d0/ 
      data c2/2.88231761d0,2.88231639d0,2.86675783d0,2.71379764d0/ 
      data c3/1.27161903d0,1.27159453d0,1.30271165d0,1.42507981d0/ 
      data c4/3.36086579d-1,3.35841662d-1,3.04724541d-1,2.55777281d-1/ 
      data c5/5.43692706d-2,5.31446837d-2,6.87032441d-2,7.84926959d-2/ 
      data c6/4.46684068d-3,2.01766691d-3,-1.09404515d-3,-1.87720131d-3/ 
      data d1/1.11841602d1,1.11841602d1,1.11823451d1,1.10708116d1/ 
      data d2/1.07900220d1,1.07900219d1,1.08009129d1,1.10685932d1/ 
      data d3/5.04405582d0,5.04405368d0,5.01682620d0,4.74914588d0/ 
      data d4/1.48355553d0,1.48352696d0,1.51983026d0,1.66259311d0/ 
      data d5/2.94075757d-1,2.93861454d-1,2.66633974d-1,2.23805121d-1/ 
      data d6/3.80584894d-2,3.72012786d-2,4.80922708d-2,5.49448872d-2/ 
      data d/2.60565706d-3,1.17697237d-3,-6.38193005d-4,-1.09503410d-3/ 
      data a/-3.53553400d-1,1.92450000d-1,-1.25000000d-1,8.94427000d-2 
     1      ,-1.76777000d-1,6.41500000d-2,-3.12500000d-2,1.78885000d-2 
     2      ,-8.83883000d-2,2.13833000d-2,-7.81250000d-3,3.57771000d-3 
     3      ,1.16317000d1,-4.41942000d-2,7.12778000d-3,-1.95313000d-3 
     4      ,7.15541750d-4/ 
      data b/6.666667d-1,8.22467d-1,7.102746d-1,6.467679d0 
     7      ,4.00000000d-1,2.46740000d0,-7.10275000d-1,-2.77186000d0 
     8      ,2.85714000d-1,4.11234000d0,3.55137000d0,2.77186000d0 
     9      ,2.22222200d-1,5.75727000d0,2.48596000d1,-6.4676800d0 
     a      ,7.79429075d-3,-4.94746507d-2,1.94857269d-2,3.56600189d-2 
     b      ,-1.73161277d-1,3.41000220d-2/ 
      data c/-7.07106800d-1,5.77350000d-1,-5.00000000d-1,4.47213500d-1 
     1      ,1.00000015d0,-4.11233500d-1,-1.77568650d0,-2.91045555d1/ 
      data ck/8.86227d-1,1.32934d0,3.32335d0/ 
      data uio/0.43139881d0,1.7597537d0,4.1044654d0,7.7467038d0, 
     u         1.3457678d1/ 
      data ui1/0.81763176d0,2.4723339d0,5.1160061d0,9.0441465d0, 
     1         1.5049882d1/ 
      data ui2/1.2558461d0,3.2070406d0,6.1239082d0,1.0316126d1, 
     2         1.6597079d1/ 
      data cio/0.37045057d0,0.41258437d0,9.777982d-2,5.3734153d-3, 
     c         3.8746281d-5/ 
      data ci1/0.39603109d0,0.69468795d0,0.2232276d0,1.5262934d-2, 
     1         1.3081939d-4/ 
      data ci2/0.76934619d0,1.7891437d0,0.70754974d0,5.6755672d-2, 
     2         5.557148d-4/ 
      data aio/0.64959979d0,0.17208724d0,0.016498837d0,4.321647d-4, 
     a         1.4302261d-6/ 
      data ai1/0.44147594d0,8.4387677d-2,5.9999383d-3,1.180802d-4, 
     1         2.9101763d-7/ 
      data ai2/0.28483475d0,4.0476222d-2,2.1898807d-3,3.3095078d-5, 
     2         6.194128d-8/ 
      data xxi/7.265351d-2,0.2694608d0,0.5331220d0,0.7868801d0, 
     x         0.9569313d0/ 
      data aai/3.818735d-2,0.1256732d0,0.1986308d0,0.1976334d0, 
     a         0.1065420d0/ 
      data cci/0.26356032d0,1.4134031d0,3.5964258d0,7.0858100d0, 
     c         1.2640801d1/ 
      data bbi/0.29505869d0,0.32064856d0,7.391557d-2,3.6087389d-3, 
     b         2.3369894d-5/ 
      data pc1/0.5d0/,pc2/0.7d0/ 
      data cpp/5.d0,1.d1,1.5d1,2.5d1,2.5d0/ 
      data fgs/0.571428571d0,0.33333333d0,0.2272727d0,0.168269d0, 
     * 0.142857143d0,5.5555555d-2,2.840909d-2,1.68269d-2/ 
      data (abc(k),k=1,76)/7.72885519d1, 
     1  1.42792768d2, 4.30552910d1, 2.43440537d0, 
     2  1.75674547d-2, 9.99400362d1, 2.73430265d2, 1.00130386d2, 
     3  6.91871969d0, 5.93132645d-2, 2.30460043d2, 7.56122303d2, 
     4  3.19543255d2, 2.57313963d1, 2.51960145d-1,-2.35425685d1, 
     5 -4.38697759d1,-1.32985534d1,-7.52438243d-1,-5.42994019d-3, 
     6 -3.05287674d1,-8.42357074d1,-3.09412790d1,-2.13850223d0, 
     7 -1.83331906d-2,-7.06062732d1,-2.33317365d2,-9.87584116d1, 
     8 -7.95332909d0,-7.78785901d-2, 1.42401164d-1, 4.12778210d-1, 
     9  1.52786427d-1, 8.84665279d-3, 6.38815164d-5, 2.18702630d-1, 
     a  8.82651141d-1, 3.60865258d-1, 2.51545288d-2, 2.15684504d-4, 
     b  5.87304073d-1, 2.59226969d0, 1.15817403d0, 9.35640728d-2, 
     c  9.16218624d-4, 2.94914091d-1, 5.29893251d-1, 1.56942521d-1, 
     d  8.85295620d-3, 6.38816670d-5, 3.77889882d-1, 1.00545595d0, 
     e  3.64435019d-1, 2.51594259d-2, 2.15684607d-4, 8.63109766d-1, 
     f  2.76526224d0, 1.16235562d0, 9.35691781d-2, 9.16218718d-4, 
     g  7.22774549d-1, 6.91700407d-1, 6.36940508d-1, 5.69038300d-1, 
     h  5.14846835d-1, 9.63560320d-1, 2.11340310d0, 4.29642580d0, 
     i  7.78581000d0, 1.33408010d1, 5.08574570d-2, 1.88622560d-1, 
     k  3.73185400d-1, 5.50816070d-1, 6.69851910d-1, 2.89632880d-1/
        data (abc(k),k=77,85)/ 
     *  4.66144392d-1, 1.53210873d-1, 1.00694874d-2, 8.53586810d-5, 
     *  1.46896384d-2, 4.60107809d-2, 6.75861819d-2, 6.21820743d-2, 
     *  3.16690579d-2/ 
          data nitm/60/ 
          data pi2/9.8696044011d0/ 
          data t5/1.3d1/,t4/1.5d1/,ro3/3.d1/,ro4/3.9d1/ 
          data nfil/1/ 
c *** 
c..        if(nfil.ne.1) go to 4000 
c..        open(unit=101,file='epeos.mes') 
c..        write(6,5001) 
c..        nfil=0 
c4000  if((nz.ge.0).and.(nz.le.5)) go to 4004 
c..      write(6,5002) nz 
c..      print 4005,nz 
c..      print 4006 
c..      stop 
4005  format('  illegal value of nz:  nz=',i5) 
4006  format(' allowed values are nz=0,1,2,3,4,5  *stop in epeos*') 
5001    format(10x,'epeos  ***  version 1.1 august 2, 1992  ***  epeos') 
5002    format(20x,'epeos  ***  error in   nz   ***  epeos'/ 
     1 1x,'illegal value of nz:  nz=',i5, 
     2 1x,'! allowed values are nz=0,1,2,3,4,5 *stop*') 
5012    format(20x,'epeos  ***  error in   t    ***  epeos'/ 
     1 1x,1p,'temperature t must be positive: t=',d13.6, 
     2 1x,'jkk=',i4,' *stop*') 
5013    format('  temperature t must be positive: t=', 
     1 1p,d13.6,' jkk=',i4,' *stop in epeos*') 
5022    format(20x,'epeos  ***     warning!     ***  epeos'/ 
     1 1x,1p,'negative or zero density: den=',d13.6,' jkk=',i4/ 
     2 1x,'calculations are going on with den=0.') 
5023    format('  negative or zero density den=', 
     1  d13.6,' jkk=',i4,' *warning in epeos* ') 
5032    format(20x,'epeos  *** error in as or zs ***  epeos'/ 
     1 1x,1p,'as and zs must be positive: as=',d10.3,' zs=',d10.3, 
     2 1x,' jkk=',i4,' *stop*') 
5033    format(' as and zs must be positive: as=',1p,d10.3, 
     1 1x,'zs=',d10.3,' jkk=',i4,'*stop in epeos*') 
c 
4004  continue 



c..use the given abar & zbar, get derivatives, entropy, include pairs
       jurs  = 0
       lst   = 2
       kentr = 1
       kpar  = 1
       eosfail = .false.

c..
c..start the pipeline loop here
      do j = jlo_eos,jhi_eos
       t   = temp_row(j) * 1.0d-9
       den = den_row(j) * 1.0d-7
       as  = abar_row(j)
       zs  = zbar_row(j)
       scn = 2.5d0 * log(abar_row(j))

c..flag set to find the working region
       nz  = 0

c..check the input 
      if(t.gt.0.d0) go to 4010 
       write(6,5012) t,jkk 
       stop 
4010  continue 
      if(pl.gt.0.d0) go to 102 
      write(6,5022) pl,jkk 
      eosfail = .true.
      stop

c      pt=3.025884d-2*t**3/cu(17) 
c      p=0.25d0*t*pt 
c      ppl=0.d0 
c      gam=cu(36) 
c      da=0.25d0 
c      go to 9 


  102 if(jurs.eq.1) stop 'tried a call to chemic'
c..call chemic 
c *** subroutine chemic calculates as,zs,scn 
      if((as.gt.0.d0).and.(zs.gt.0.d0)) go to 4030 
      write(6,5032) as,zs,jkk 
      eosfail = .true. 
      stop 
c..
c..start normal calculations
4030  emue=as/zs 
      rg=cu(2)/emue 
      ki=1 
      hpr=0.d0 
   90 alf=cu(1)/t 
      al1=cu(4)/alf 
      plm=pl/emue 
      sqe=0.d0 
      ei=t*rg 
      pi=ei*pl 
      eg=cu(3)*ei 
  590  continue 
c 
c *** search for required working region 
      if(nz.ne.0) go to(1,2,3,4,5),nz 
      if(ki.ne.1) go to  123 
      if(plm.le.ro3) go to 310 
      if(plm.ge.ro4) go to 4 
      if(t.le.t5) go to 550 
      if(t.ge.t4) go to 4 
c *** searching around the triangle 
          x=(ro4-ro3)/(t4-t5) 
          y=ro3+(t4-t)*x 
          if(y.gt.plm) go to 800 
          go to 4 
  310  if(t.le.t5) go to 123 
          if(t.ge.t4) go to 4 
c *** interpolation over t in region 45 
          t1=t5 
          t2=t4 
          nz2=4 
          nz1=5 
          go to 128 
  550  continue 
c *** interpolation over density for density < ro4 
      nzp1=0 
      if(t.lt.tt2) nzp1=3 
  583  kin=ki 
       nz=nzp1 
       ki=6 
       go to 590 

  577  pn1=pe 
       en1=ee 
       sn1=se 
       psn1=psi 
       hprn=hpr 
       pnp1=ppl 
       enp1=epl 
       snp1=spl 
       pnt=pt 
       ent=et 
       snt=st 
       nzr1=nzr 
       nz=4 
       ki=7 
       go to 590 
  578  z1=ro4-ro3 
       x2=(plm-ro3)/z1 
       x=x2**2*(3.d0-2.d0*x2) 
       x1=1.d0-x 
       z1=z1*emue 
       x3=6.d0*x2*(1.d0-x2)/z1 
       p=pe 
       e=ee 
       s=se 
       pe=pe*x+pn1*x1 
       ee=ee*x+en1*x1 
       if(kentr.eq.0) go to 591 
       se=se*x+sn1*x1 
  591 if(lst.eq.0) go to 592 
      pt=pt*x+pnt*x1 
      et=et*x+ent*x1 
      ppl=ppl*x+pnp1*x1+(p-pn1)*x3 
      epl=epl*x+enp1*x1+(e-en1)*x3 
      if(kentr.eq.0) go to 592 
      st=st*x+snt*x1 
      spl=spl*x+snp1*x1+(s-sn1)*x3 
  592  psi=psi*x+psn1*x1 
      hpr=hpr*x+hprn*x1 
      nzr=10*nzr1+nzr 
      ki=kin 
      go to 134 
  800 continue 
c ***** the triangle 
          nz2=4 
         nz1=0 
          t1=t5 
          t2=t4-(plm-ro3)/x 
          go to 128 
  123 if(ki.ne.4) go to 136 
      if(nz2.eq.5) go to 111 
      nzp1=5 
          go to 583 
  136 if(plm.lt.ps1) go to 110 
      if(t.lt.tt1) go to 111 
      kzz=2 
         go to 121 
  115 y=x*grcf 
      if(plm.gt.y) go to 3 
      kzz=3 
         z1=t 
      if(t.lt.tt2) go to 112 
      if(plm.lt.x) go to 5 
  113 kkz=0 
  116 go to 121 
  114 dl=(x-plm)/xz 
      x1=1.d0 
      if(dl.lt.0.d0) x1=-1.d0 
      x2=dl*x1 
      if(x2.gt.0.3d0) dl=0.3d0*x1 
      t=t*(1.d0-dl) 
      if(x2.gt.eit)go to 116 
      if(kkz.eq.1) go to 118 
      t2=t 
      t=z1 
  138 z2=plm 
         plm=plm/grcf 
         kkz=1 
      go to 116 
  118 t1=t 
         nz1=3 
         plm=z2 
         t=z1 
         nz2=5 
      go to 128 
  112 if(plm.gt.pss1) go to 117 
      t1=tt1 
         t2=tt2 
        nz1=0 
        nz2=5 
      go to 128 
  117 if(plm.gt.pss2) go to 113 
      t2=tt2 
      go to 138 
  110 if(t.lt.tpar) go to 111 
      if(t.gt.tt2) go to 5 
      sqe=exp(-alf) 
      if(t.gt.tt1) go to 119 
      if(plm.gt.ro1*sqe) go to 1 
      if(t.gt.tpar1) go to 119 
      if(plm.gt.rotp) go to 120 
      t1=tpar 
      t2=tpar1 
  122 nz1=1 
      nz2=5 
      sqe=0.d0 
      go to 128 
  119 if(plm.lt.ro2*sqe) go to 5 
  120 t1=log(plm) 
      t2=cu(1)/(ro2l-t1) 
      t1=cu(1)/(ro1l-t1) 
      go to 122 
  111 sq=t*sqrt(t) 
      y=2.095d-2*sq 
      x=y*grif 
      if(plm-x)44,44,51 
c 
c perfect gas with corrections for degeneracy and pairs (nz=1) 
    1 sq=t*sqrt(t) 
   44 nzr=1 
      qa=6.9712909d0*plm/sq 
      x=qa*al1 
      x2=qa-x*(cu(5)-cu(6)*al1) 
      pe=cu(4)+x2 
      epl=qa-x*(cu(10)+cu(11)*al1) 
      x8=cu(9)*al1 
      x3=x8*(cu(4)+(cu(14)*al1-cu(4))*al1) 
      ee=cu(4)+epl+x3 
      if(lst.eq.0) go to 6 
      pt=cu(4)-cu(16)*x2-x*(cu(5)-cu(15)*cu(6)*al1) 
      ppl=cu(4)+cu(15)*x2 
      et=cu(4)-cu(16)*epl+x8*(cu(15)+al1*(cu(18)*al1-cu(17))-qa* 
     1   (cu(19)+cu(20)*al1)) 
    6 if(t.lt.tpar) go to 8 
      if(sqe.eq.0.d0) sqe=exp(-alf) 
      x4=0.268192d0*(al1*sqe/plm)**2 
      x5=x4*al1 
      hpr=x5*(cu(4)+al1*(cu(7)+cu(8)*al1)) 
      pe=pe+hpr 
      x6=(x4+x5*(cu(12)+cu(13)*al1))/cu(3) 
      ee=ee+x6 
      if(lst.eq.0) go to 8 
      ppl=ppl-hpr 
      pt=pt+hpr*cu(15)*(cu(15)+alf)+x5*(cu(7)+cu(21)*al1)*al1 
      x8=x6*cu(15) 
      et=et+x8*(cu(3)+alf)+x5*(cu(23)+cu(22)*al1) 
      epl=epl-x8 
    8 if((kentr+ki).eq.1) go to 56 
      x7=cu(16)*(qa+x*(cu(5)-cu(25)*al1)) 
      psi=log(cu(29)*qa) 
      se=cu(24)-psi+x7+al1*(cu(7)+al1*(al1*cu(26)-cu(5))) 
      psi=psi+2.d0*x2+0.5d0*hpr-1.875d0*al1* 
     * (1.d0+al1*(al1*0.1875d0-0.5d0)) 
      if(kentr.eq.0) go to 56 
      if(lst.eq.0) go to 53 
      st=cu(3)*(cu(4)+al1*(cu(27)+al1*(cu(5)*al1-cu(7)))) 
      st=st-cu(28)*qa*(cu(17)+al1*(cu(5)+cu(25)*al1)) 
      spl=x7-cu(4) 
   53 if(t.lt.tpar) go to 101 
      x8=x4+x5*(cu(30)+cu(31)*al1) 
      se=se+x8 
      if(lst.eq.0) go to 10 
      spl=spl-cu(15)*x8 
      st=st+x4*(cu(15)*alf+cu(34)+al1*(cu(32)+cu(33)*al1)) 
  101 if(lst.eq.0) go to 10 
      st=st*rg/t 
      spl=spl*rg 
   10 se=se*rg 
   56 pe=pi*pe 
      ee=eg*ee 
      if(lst.eq.0) go to 50 
      pt=rg*pl*pt 
      et=eg*et/t 
      ppl=ei*ppl 
      epl=eg*epl 
   50 go to 135 
c 
c *** addition the ion and black-body radiation components to eos 
   57 continue 
c ********************************************************************* 
      x=pi/zs 
      x1=eg/zs 
      v=7.56471d-3*t**4 
c      v = 7.56590624067000303d-3 * t**4
      x3=v/cu(17) 
      p=x+pe+x3 

      pion = x
      dpiondd = x/pl
      dpiondt = x/t


      prad = x3


      beta=x3/p 
      x3=v/pl 
      e=x1+x3+ee
 
      eion = x1
      erad = x3 

      if(kentr.eq.0) go to 7 
      x6=cu(2)/as 
      x4=pl*(cu(35)/(t*sqrt(t))) 
      x4=log(x4) 
      x4=x6*(cu(24)-x4+scn) 
      x5=cu(36)*x3/t 
      s=x5+x4+se 

      sion = x4
      srad = x5

      sek=se/cu(2) 
      sk=s/cu(2) 
      if(lst.eq.0) go to 9 

c store the electron-positron entropy derivatives before adding rad+ion
      dsepdd = spl
      dsepdt = st
      st=st+(cu(17)*x5+x1/t)/t 
      spl=spl-x6-x5 

c..dsraddt and dsraddd 
c      write(6,*) cu(17)*x5/t*1.0d-1,-x5*1.0d1/den

c..dsiondt and dsiondd
c      write(6,*) x1/t/t * 1.0d-1,-x6*1.0d1/den


      go to 45 
    7 if(lst.eq.0) go to 9 
   45 continue

c..store the electron-positron pressure derivatives before adding rad+ion
      dpepdd = ppl
      dpepdt = pt
      pt=pt+(x+cu(36)*v)/t 
      ppl=ppl+x/pl 

c..store the electron-positron energy derivatives before adding rad+ion
      deepdd = epl
      deepdt = et
      et=et+(x1+cu(37)*x3)/t 
      epl=epl-x3 

c ********************************************************************* 
      x4=pt/(pl*et) 
      x5=pl/p 
      gam=x5*(ppl+t*x4*pt/pl) 
      da=x4/gam 
      cp=gam*et*(p/ppl) 
      dpe=x5*(epl+t*pt/pl)-1.d0 
      if(kentr.ne.0) then 
      dse=t*st/et-1.d0 
      dsp=-spl*(pl/pt)-1.d0 
      endif 


c..normal exit from epeos 
c..set the return arguments in cgs units
   9  continue

c  *** dpe = (den/p)(epl+t*pt/den)-1 -- thermodynamic identity: dpe=0 
c  *** dse = t*st/et-1 -- thermodynamic identity: dse=0 
c  *** dsp = -spl*(den/pt)-1 -- thermodynamic identity: dsp=0 
c       dse = t*st/et - 1.0d0
c       dpe = (epl*den + t*pt)/p - 1.0d0
c       dsp = -spl*(den/pt) -1.0d0
c       write(6,696) dse,dpe,dsp

      p       = p       * 1.0d24
      prad    = prad    * 1.0d24
      pion    = pion    * 1.0d24
      pe      = pe      * 1.0d24

      pt      = pt      * 1.0d15
      dpiondt = dpiondt * 1.0d15
      dpepdt  = dpepdt  * 1.0d15

      ppl     = ppl     * 1.0d17
      dpiondd = dpiondd * 1.0d17
      dpepdd  = dpepdd  * 1.0d17

      e       = e       * 1.0d17
      erad    = erad    * 1.0d17
      eion    = eion    * 1.0d17
      ee      = ee      * 1.0d17

      et      = et      * 1.0d8
      deepdt  = deepdt  * 1.0d8

      epl     = epl     * 1.0d10/den
      deepdd  = deepdd  * 1.0d10/den

      s       = s       * 1.0d8
      sion    = sion    * 1.0d8 
      srad    = srad    * 1.0d8 
      se      = se      * 1.0d8

      st      = st      * 1.0d-1
      dsepdt  = dsepdt  * 1.0d-1

      spl     = spl     * 1.0d1/den
      dsepdd  = dsepdd  * 1.0d1/den


c..coulomb section:
c..initialize

       pcoul    = 0.0d0
       dpcouldd = 0.0d0
       dpcouldt = 0.0d0
       ecoul    = 0.0d0
       decouldd = 0.0d0
       decouldt = 0.0d0
       scoul    = 0.0d0
       dscouldd = 0.0d0
       dscouldt = 0.0d0


c..uniform background corrections only 
c..see yakovlev & shalybkov 1989 
c..lami is the average ion seperation
c..plasg is the plasma coupling parameter
       z        = forth * pie
       xni      = avo/as * den_row(j)
       dxnidd   = avo/as
       sss      = z * xni
       dsdd     = z * dxnidd
       lami     = 1.0d0/sss**third
       inv_lami = 1.0d0/lami
       z        = -third * lami/sss 
       lamidd   = z * dsdd
       ktinv    = 1.0d0/(kerg*temp_row(j))
       plasg    = zs*zs*esqu*ktinv*inv_lami
       z        = -plasg * inv_lami 
       plasgdd  = z * lamidd
       plasgdt  = -plasg*ktinv * kerg

c..fermi temperature, melting temperature and density condition
       y        = den_row(j) * 1.0d-6
       z        = 1.009d0 * (y * zs/as)**third
       tfermi   = 5.930d9 * (sqrt(1.0d0 + z*z) - 1.0d0)
       tmelt    = 1.278d5 * zs*zs*(y/as)**third
       rhocond  = as*zs


c..for any of these uniform background corrections to apply
c..the following conditions must hold
c..temp << t_fermi   electron gas is degenerate
c..temp >  t_melt    temperature exceeds the ion melting temperature
c..rho >>  abar*zbar electron gas is almost perfect, ionization full


       if ( (temp_row(j) .lt. tfermi) .and.
     1      (temp_row(j) .gt. 0.1d0 * tmelt)  .and.
     2      (den_row(j)  .gt. 2.0d0 * rhocond)     
     3       ) then


c..yakovlev & shalybkov 1989 equations 82, 85, 86, 87
        if (plasg .ge. 1.0) then
         x        = plasg**(0.25d0) 
         y        = avo/as * kerg 
         ecoul    = y * temp_row(j) * (aa1*plasg + bb1*x + cc1/x + dd1)
         pcoul    = third * den_row(j) * ecoul
         scoul    = -y * (3.0d0*bb1*x - 5.0d0*cc1/x
     1              + dd1 * (log(plasg) - 1.0d0) - ee1)

         y        = avo/as*kerg*temp_row(j)
     1              * (aa1 + 0.25d0/plasg*(bb1*x - cc1/x))
         decouldd = y * plasgdd 
         decouldt = y * plasgdt + ecoul/temp_row(j)

         y        = third * den_row(j)
         dpcouldd = third * ecoul + y*decouldd
         dpcouldt = y * decouldt

         y        = -avo*kerg/(as*plasg)
     1                *(0.75d0*bb1*x +1.25d0*cc1/x +dd1)
         dscouldd = y * plasgdd
         dscouldt = y * plasgdt


c..yakovlev & shalybkov 1989 equations 102, 103, 104
        else if (plasg .lt. 1.0) then
         x        = plasg*sqrt(plasg)
         y        = plasg**bb2
         z        = cc2 * x - third * aa2 * y
         pcoul    = -pion * z
         ecoul    = 3.0d0 * pcoul/den_row(j)
         scoul    = -avo/as*kerg*(cc2*x -aa2*(bb2-1.0d0)/bb2*y)

         sss      = 1.5d0*cc2*x/plasg - third*aa2*bb2*y/plasg
         dpcouldd = -dpiondd*z - pion * sss * plasgdd
         dpcouldt = -dpiondt*z - pion * sss * plasgdt

         sss      = 3.0d0/den_row(j)
         decouldd = sss * dpcouldd - ecoul/den_row(j)
         decouldt = sss * dpcouldt

         sss      = -avo*kerg/(as*plasg)
     1               * (1.5d0*cc2*x -aa2*(bb2-1.0d0)*y)
         dscouldd = sss * plasgdd
         dscouldt = sss * plasgdt
        end if


c..soften the edges of the coulomb correction turn-on
c..by doing linear interpolation to zero coulomb corrections
c..within reasonable limits of when the above expressions are valid

        z2 = 1.0d0
        if (den_row(j) .gt. 2.0d0*rhocond .and. 
     1      den_row(j) .lt. 20.0d0*rhocond) then
         x        = den_row(j)
         x1       = 2.0d0 * rhocond
         x2       = 20.0d0 * rhocond
         z2     = (x-x1)/(x2-x1)

        else if (temp_row(j) .gt. 0.1d0 * tmelt  .and. 
     1           temp_row(j) .lt. tmelt) then
         x        = temp_row(j)
         x1       = 0.1d0 * tmelt
         x2       = tmelt
         z2     = (x-x1)/(x2-x1)

        else if (temp_row(j) .gt. 0.5d0 * tfermi  .and. 
     1           temp_row(j) .lt. tfermi) then
         x        = temp_row(j)
         x1       = tfermi
         x2       = 0.5d0 * tfermi
         z2     = (x-x1)/(x2-x1)
        end if

        pcoul    = z2 * pcoul
        dpcouldd = z2 * dpcouldd
        dpcouldt = z2 * dpcouldt
 
        ecoul    = z2 * ecoul
        decouldd = z2 * decouldd
        decouldt = z2 * decouldt
 
        scoul  =   z2 * scoul
        dscouldd = z2 * dscouldd
        dscouldt = z2 * dscouldt
       end if





c..and add in the coulomb corrections
      p   = p + pcoul
      pt  = pt + dpcouldt
      ppl = ppl + dpcouldd

      e   = e + ecoul
      et  = et + decouldt
      epl = epl + decouldd

      sk  = sk + scoul/(kerg*avo)
      s   = s  + scoul
      st  = st + dscouldt
      spl = spl + dscouldd



c..the temperature and density exponents (c&g 9.81 9.82) 
c..the specific heat at constant volume (c&g 9.92)
c..the third adiabatic exponent (c&g 9.93)
c..the first adiabatic exponent (c&g 9.97) 
c..the second adiabatic exponent (c&g 9.105)
c..the specific heat at constant pressure (c&g 9.98) 
c..and relativistic formula for the sound speed (c&g 14.29)
      zz    = p/den_row(j)
      chit  = temp_row(j)/p * pt
      chid  = ppl/zz
      cv    = et
      x     = zz * chit/(temp_row(j) * cv)
      gam3  = x + 1.0d0
      gam1  = chit*x + chid
      nabad = x/gam1
      gam2  = 1.0d0/(1.0d0 - nabad)
      cp    = cv * gam1/chid
      z     = 1.0d0 + (e + con2)/zz
      sound = clight * sqrt(gam1/z)

c..these are the 3 thermodynamic consistency checks
c..each is zero if the consistency is perfect
        dse = temp_row(j)*st/et - 1.0d0
        dpe = (epl*den_row(j)*den_row(j) + temp_row(j)*pt)/p - 1.0d0
        dsp = -spl*(den_row(j)**2/pt) - 1.0d0

c..store this row and go to the end of the do loop
        ptot_row(j)   = p
        dpt_row(j)    = pt
        dpd_row(j)    = ppl
        dpa_row(j)    = 0.0d0
        dpz_row(j)    = 0.0d0

        etot_row(j)   = e
        det_row(j)    = et
        ded_row(j)    = epl
        dea_row(j)    = 0.0d0
        dez_row(j)    = 0.0d0

        stot_row(j)   = s
        dst_row(j)    = st
        dsd_row(j)    = spl
        dsa_row(j)    = 0.0d0
        dsz_row(j)    = 0.0d0

        prad_row(j)   = prad
        erad_row(j)   = erad
        srad_row(j)   = srad

        pion_row(j)   = pion
        eion_row(j)   = eion
        sion_row(j)   = sion
        xni_row(j)    = avo/as * den_row(j)

        pele_row(j)   = pe
        ppos_row(j)   = 0.0d0
        dpept_row(j)  = dpepdt
        dpepd_row(j)  = dpepdd
        dpepa_row(j)  = 0.0d0
        dpepz_row(j)  = 0.0d0

        eele_row(j)   = ee
        epos_row(j)   = 0.0d0
        deept_row(j)  = deepdt
        deepd_row(j)  = deepdd
        deepa_row(j)  = 0.0d0
        deepz_row(j)  = 0.0d0

        sele_row(j)   = se
        spos_row(j)   = 0.0d0
        dsept_row(j)  = dsepdt
        dsepd_row(j)  = dsepdd
        dsepa_row(j)  = 0.0d0
        dsepz_row(j)  = 0.0d0

        xnem_row(j)   = zs * xni_row(j)
        xne_row(j)    = xnem_row(j) + 0.5d0 * hpr * xnem_row(j)
        dxnet_row(j)  = 0.0d0
        dxned_row(j)  = avo/emue * (1.0d0 + hpr)
        dxnea_row(j)  = 0.0d0
        dxnez_row(j)  = 0.0d0
        xnp_row(j)    = 0.5d0 * hpr * xnem_row(j)

        etaele_row(j) = psi
        detat_row(j)  = 0.0d0
        detad_row(j)  = 0.0d0
        detaa_row(j)  = 0.0d0
        detaz_row(j)  = 0.0d0
        etapos_row(j) = 0.0d0
        beta          = kerg * temp_row(j)/mecc
        if (beta .gt. 0.02) etapos_row(j) = -psi - 2.0d0/beta

        pcou_row(j)   = pcoul
        ecou_row(j)   = ecoul
        scou_row(j)   = scoul
        plasg_row(j)  = plasg

        dse_row(j)    = dse
        dpe_row(j)    = dpe
        dsp_row(j)    = dsp

        cv_row(j)     = cv
        cp_row(j)     = cp
        gam1_row(j)   = gam1
        gam2_row(j)   = gam2
        gam3_row(j)   = gam3
        cs_row(j)     = sound

      goto 6969

c 
c *** further search for working region 
   51 if(plm.lt.y) go to 52 
      kzz=1 
  121 z=t/ep2 
      if(z.gt.0.57142857d0) go to 64 
      x=xep2*sq 
      if(kzz.eq.3) xz=cu(3)*x 
      go to 125 
   64 if(z.gt.1.07142857d0) go to 54 
      xz=7.306392d-2*z 
      x=xz+3.414385d-2 
      go to 125 
   54 if(z.gt.1.d1) go to 58 
      x=((4.77856d-2*z-1.41839d-2)*z+7.2001d-2)*z-7.20897d-3 
      if(kzz.eq.3) xz=((1.433568d-1*z-2.83678d-2)*z+7.2001d-2)*z 
      go to 125 
   58 x=4.708d-2*z**3 
      if(kzz.eq.3) xz=cu(17)*x 
  125 go to (55,115,114),kzz 
   55 if(plm-x) 59,59,61 
c 
c *** expansion over half-integer fermi--dirac functions (nz=2) 
    2 sq=t*sqrt(t) 
   59 nzr=2 
      kenf=1 
      zp=plm/(cu(38)*sq) 
      x=psi 
      nit=0 
      kf=1 
      x1=cu(9)*al1 
      al2=al1**2 
      al3=0.21875d0*al2 
      go to 26 
   11 z=zp-f12-x1*f32-al3*f52 
      v=f12s+x1*f32s+al3*f52s 
      dl=z/v 
      z1=x 
      if(z1.lt.0.d0)z1=-x 
      v=1.d0 
      if(dl.lt.0.d0) v=-1.d0 
      z=v*dl 
      if(z1.lt.0.1d0) go to 41 
      z=z/z1 
      if(z.gt.0.3d0) dl=0.3d0*v*z1 
   41 x=x+dl 
      nit=nit+1 
      if(z.lt.eit) go to 73 
      if(nit.lt.nitm) go to 26 
   66  write(6,5072) nzr 
       write(6,5073) dl,x,eit,nitm,jkk 
c..       print 72,nzr 
c..       print 4072,dl,x,eit,nitm,jkk 
       eosfail = .true. 
       stop 
   72  format('  iterations in region nzr=',i1, 
     1 1x,'do not converge!') 
4072   format('  dx=',1p,d10.3,' x=',d10.3,' eit=',d10.3, 
     1 1x,'nitm=',i3,' jkk=',i4,' *stop in epeos*') 
5072   format(10x,'epeos  ***  runaway of iterations in region nzr=', 
     1 i1,' *** stop*') 
5073   format(10x,1p,'dx=',d10.3,' x=',d10.3,' eit=',d10.3,' nitm=',i3, 
     1 1x,'jkk=',i4) 
c 
   73 kf=2 
      go to 26 
  300 go to(31,9),kenf 
   31 psi=x 
      v=sq*al1 
      x=cu(19)*al1 
      z=cu(39)*v 
      z1=cu(3)*z/pl 
      z3=x1*f32 
      z4=x*f52 
      al4=9.375d-2*al2 
      y=al4*f72 
      pe=z*(f32+z4+y) 
      z5=x1*f52 
      y1=al3*f72 
      ee=z1*(f32+z5+y1) 
      if(kentr.eq.0) go to 34 
      z2=cu(41)*sq/pl 
      dl=cu(40)*psi 
      z10=cu(44)*psi 
      z11=cu(45)*al1 
      z6=z11*(f52-dl*f32) 
      al5=0.16875d0*al2 
      y3=0.77777777778d0*psi 
      y2=al5*(f72-y3*f52) 
      se=z2*(f32-z10*f12+z6+y2) 
   34 if(lst.eq.0) go to 35 
      pal=f12s+x1*f32s+al3*f52s 
      pap=zp/pal 
      pal=(cu(3)*zp+z3+cu(15)*al3*f52)/pal 
      z8=f32s+x1*f52s+al3*f72s 
      z7=z5+cu(15)*y1-pal*z8 
      et=(cu(27)*ee+z1*z7)/t 
      epl=z1*z8*pap-ee 
      z9=f32s+x*f52s+al4*f72s 
      pt=(cu(27)*pe+z*(z4+cu(15)*y-pal*z9))/t 
      ppl=z*z9*pap/pl 
      if(kentr.eq.0) go to 35 
      z9=f32s-z10*f12s-cu(44)*f12+z11*(f52s-dl*f32s-cu(40)*f32) 
      z9=z9+al5*(f72s-0.77777777778d0*f52-y3*f52s) 
      st=(cu(3)*se+z2*(z6+cu(15)*y2-pal*z9))/t 
      spl=z2*pap*z9-se 
   35 go to 135 
c 
c *** procedure of calculation of half-integer fermi--dirac functions 
c..      entry fd12f 
c      kenf=2 
c      kf=2 
c..      x=psi 
   26 if(x.ge.-1.d0) go to 21 
      z=exp(x) 
      v=ck(1)*z 
      f12=v*(1.d0+z*(a(1)+z*(a(2)+z*(a(3)+z*a(4))))) 
      f12s=v*(1.d0+z*(c(1)+z*(c(2)+z*(c(3)+z*c(4))))) 
      f32=ck(2)*z*(1.d0+z*(a(5)+z*(a(6)+z*(a(7)+z*a(8))))) 
      f52=ck(3)*z*(1.d0+z*(a(9)+z*(a(10)+z*(a(11)+z*a(12))))) 
      go to(30,12),kf 
   12 f72=a(13)*z*(1.d0+z*(a(14)+z*(a(15)+z*(a(16)+z*a(17))))) 
      go to 30 
   21 if(x.ge.-.1d0) go to 22 
      n=1 
      go to 14 
   22 if(x.ge.1.d0) go to 23 
      n=2 
      go to 14 
   23 if(x.ge.2.5d0) go to 24 
      n=3 
      go to 14 
   24 if(x.ge.4.5d0) go to 25 
      n=4 
   14 f12=a1(n)+x*(a2(n)+x*(a3(n)+x*a4(n))) 
      f12s=a2(n)+x*(df1(n)+x*df2(n)) 
      f32=b1(n)+x*(b2(n)+x*(b3(n)+x*(b4(n)+x*b5(n)))) 
      f52=c1(n)+x*(c2(n)+x*(c3(n)+x*(c4(n)+x*(c5(n)+x*c6(n))))) 
      go to(30,13),kf 
   13 f72=d1(n)+x*(d2(n)+x*(d3(n)+x*(d4(n)+x*(d5(n)+x*(d6(n)+x*d(n)))))) 
   30 f32s=1.5d0*f12 
      f52s=2.5d0*f32 
      if(kf.eq.1) go to 11 
      f72s=3.5d0*f52 
      go to 300 
   25 z=sqrt(x) 
      z1=z*x 
      z2=1.d0/x**2 
      f12=z1*(b(1)+(b(2)+(b(3)+b(4)*z2)*z2)*z2) 
      f12s=z*(c(5)+(c(6)+(c(7)+c(8)*z2)*z2)*z2) 
      z=z1*x 
      f32=z*(b(5)+(b(6)+(b(7)+b(8)*z2)*z2)*z2) 
      f52=x*z*(b(9)+(b(10)+(b(11)+b(12)*z2)*z2)*z2) 
      f32=f32+b(17) 
      f52=f52+b(18)+b(19)*x 
      go to(30,17),kf 
   17 f72=z*(b(13)+(b(14)+(b(15)+b(16)*z2)*z2)*z2)/z2 
      f72=f72+b(20)+x*(b(21)+b(22)*x) 
      go to 30 
c 
c *** search for working region is continued 
   61 y=x*grcf 
      if(plm.lt.y) go to 62 
c 
c *** chandrasekhar's expansion for degenerate gas (nz=3) 
    3 nzr=3 
      x1=plm/cu(47) 
      if(psi.lt.3.d0) psi=3.d0 
      nit=0 
      x=psi*al1 
      x=sqrt(x*(x+cu(15))) 
      al2=al1**2 
      z2=1.644934d0*al2 
      z4=1.894066d0*al2**2 
   74 x2=x**2 
      x3=x2*x 
      x5=cu(17)*(z4/x3)/x2 
      x6=z2/x 
      x8=cu(15)*z2*x 
      dl=(x3*cu(49)+x8+x6+x5-x1)/(x3+x8-x6-cu(50)*x5) 
      x6=1.d0 
      if(dl.lt.0.d0) x6=-1.d0 
      z1=x6*dl 
      if(z1.gt.0.9d0) dl=0.9d0*x6 
      x=x*(cu(4)-dl) 
      nit=nit+1 
      if(z1.lt.eit) go to 71 
      if(nit.lt.nitm) go to 74 
      go to 66 
   71 x2=x**2 
      x4=cu(4)+x2 
      z=sqrt(x4) 
      y1=cu(4)+z 
         z1=x2/y1 
      psi=alf*z1 
      z3=x*z 
      x3=x*x2 
      z5=cu(15)*x2 
      x7=z5+cu(4) 
      if(x.gt.0.1d0) go to 174 
      x5=x2*x3 
      f0=x5*(1.6d0-x2*(fgs(1)-x2*(fgs(2)-x2*(fgs(3)-x2*fgs(4)))))*cu(49) 
      g0=x5*(0.8d0-x2*(fgs(5)-x2*(fgs(6)-x2*(fgs(7)-x2*fgs(8))))) 
      go to 175 
  174 x6=log(x+z) 
      f0=z3*(z5-cu(17))*cu(49)+x6 
      g0=x7*z3-x6-x3*cu(52) 
  175 x5=z4/x3 
      y5=z5-cu(4) 
      f2=z2*cu(51)*z3 
      f4=x5*cu(51)*y5*z 
      pe=cu(46)*(f0+f2+f4) 
      y2=cu(51)/y1 
      y4=z*y1 
      g2=z2*y2*x*(x7+y4) 
      g4=x5*cu(17)*y2*(cu(4)+y5*y4) 
      ee=cu(46)*(g0+g2+g4)/pl 
      if(kentr.eq.0) go to 75 
      y6=cu(48)/(t*pl) 
      se=y6*(f2+cu(15)*f4) 
   75 if(lst.eq.0) go to 76 
      z6=cu(54)*x5/x3 
      z7=x7-cu(15) 
      pap=x2+z2*z7/x2-z6 
      pal=cu(15)*(z2*x7+cu(55)*x5/x)/(x*pap) 
      pap=x1/pap 
      z9=cu(51)*x/z 
      z10=x1*z9 
      z11=cu(46)*pap/pl 
      ppl=z11*z10 
      y3=cu(46)/t 
      pt=y3*(cu(15)*(f2+cu(15)*f4)-z10*pal) 
      g0=z1*x2*cu(51) 
      v=cu(15)*(g2+cu(15)*g4) 
      g2=cu(51)*z2*((cu(4)+cu(17)*x7*y1)/y4-cu(15)) 
      g4=cu(53)*x5*(cu(37)-z)/(x*y4) 
      g0=g0+g2+g4 
      epl=z11*g0-ee 
      et=y3*(v-pal*g0)/pl 
      if(kentr.eq.0) go to 76 
      g4=(z2*x7+cu(55)*x5/x)*z9/x 
      spl=y6*pap*g4-se 
      st=(se+y6*(cu(37)*f4-pal*g4))/t 
   76 go to 135 
c 
c *** quadratures taken with the gauss method (nz=5) 
    5 nzr=5 
      kkk=1 
         kk1=1 
      al3=alf**3 
         z11=plm*al3 
         x2=cu(15)*alf 
      z10=z11/cu(47) 
      kw=1 
         ku=10 
      go to 151 
  169 z=(g1m-z10)/g1mp 
      z3=1.d0 
      if(z.lt.0.d0) z3=-1.d0 
      z2=z*z3 
      z1=pc 
      if(pc.lt.0.d0) z1=-pc 
      if(z1.lt.1.d0) go to 181 
      z2=z2/z1 
      if(z2.gt.0.3d0) z=0.3d0*z1*z3 
  181 pc=pc-z 
      if(z2.gt.eit) go to 151 
      ku=15 
         kw=2 
      go to 151 
  170 z=1.44059d0/(al3*alf) 
      z1=z/pl 
        z2=z/cu(17) 
      pe=z2*gp 
         ee=z1*ge 
         hpr=cu(15)*g1/z10 
      if(kentr.eq.0) go to 182 
      z3=pc+alf 
         y3=g3+g31+cu(49)*gp-z3*g1m 
      se=z1*y3/t 
  182 if(lst.eq.0) go to 183 
      pap=g1m/g1mp 
         x=g1a1-g1a 
      pal=(cu(17)*g1m-alf*(x+cu(15)*g1p))/g1mp 
      x1=g2p1-g2p 
         y2=g2a+g2a1-cu(15)*g2p 
      ppl=ei*cu(49)*x1/g1mp 
      pt=pe*(cu(37)-(alf*y2+x1*pal)/gp)/t 
      y1=g4p1-g4p-x2*g1p 
      epl=ee*(pap*y1/ge-cu(4)) 
      et=ee*(cu(37)-(alf*(g4a1+g4a)+x2*(g1-g4p+alf*g1a-x2*g1p)+y1*pal)/ 
     *ge)/t 
      if(kentr.eq.0) go to 183 
      y1=g3p1-g3p+cu(49)*x1-g1m-z3*g1mp 
      spl=se*(pap*y1/y3-cu(4)) 
      st=g3a1+g3a+y2*cu(49)-g1m-z3*(g1a1-g1a+cu(15)*g1p)-cu(15)*g3p 
      st=se*(cu(17)-(pal*y1+st*alf)/y3)/t 
  183 go to 135 
  151 kpg=0 
  152 wo=0.d0 
         w1=0.d0 
         w2=0.d0 
      wop=0.d0 
         w1p=0.d0 
        w2p=0.d0 
      woa=0.d0 
         w1a=0.d0 
        w2a=0.d0 
      if(pc.gt.pc2) go to 155 
      if(kkk.eq.0) go to 158 
      do 157 k=1,15 
  157 uac(k+65)=sqrt(uac(k)+x2) 
      kkk=0 
  158 if(pc.gt.pc1) go to 156 
      if(pc.lt.-4.4d1) go to 163 
      x=exp(pc) 
      do 161 k=1,ku 
  161 uac(k+80)=uac(k+15)/(uac(k+30)*x+1.d0) 
      do 162 k=1,5 
      z=wk1(k)*wk4(k) 
      wo=wo+z 
      wop=wop+z/(1.d0+aio(k)*x) 
      z=wk2(k)*wk5(k) 
      w1=w1+z 
      w1p=w1p+z/(1.d0+ai1(k)*x) 
      if(kw.eq.1) go to 162 
      z=wk6(k)*wk3(k) 
      w2=w2+z 
      if(lst.eq.0) go to 162 
      woa=woa+wk4(k)/wk1(k) 
      w1a=w1a+wk5(k)/wk2(k) 
      w2p=w2p+z/(1.d0+ai2(k)*x) 
      w2a=w2a+wk6(k)/wk3(k) 
  162 continue 
      wop=wop*x 
         w1p=w1p*x 
      wo=wo*x 
        w1=w1*x 
      if(kw.eq.1) go to 163 
      w2=w2*x 
      if(lst.eq.0) go to 163 
      w1a=w1a*x 
         woa=woa*x 
      w2a=w2a*x 
         w2p=w2p*x 
  163 g1=w1+alf*wo 
      g1p=w1p+alf*wop 
      if(kw.eq.1) go to 164 
      g2=w2+x2*w1 
         g3=w2+alf*(w1+g1) 
        g4=w2+alf*w1 
      if(lst.eq.0) go to 164 
      g1a=w1a+wo+alf*woa 
      g2p=w2p+x2*w1p 
         g3p=w2p+alf*(w1p+g1p) 
         g4p=w2p+alf*w1p 
      g2a=w2a+2.d0*w1+x2*w1a 
         g3a=w2a+w1+g1+alf*(w1a+g1a) 
      g4a=w2a+w1+alf*w1a 
  164 if(kpg.eq.1) go to 166 
      g11=g1 
         g1p1=g1p 
      if(kw.eq.1) go to 168 
      g21=g2 
         g31=g3 
        g41=g4 
      if(lst.eq.0) go to 168 
      g1a1=g1a 
         g2p1=g2p 
         g3p1=g3p 
         g4p1=g4p 
      g2a1=g2a 
         g3a1=g3a 
         g4a1=g4a 
  168 pc=-pc-x2 
         kpg=1 
      if(pc.gt.-4.4d1) go to 152 
         g1=0.d0 
         g2=0.d0 
         g3=0.d0 
         g4=0.d0 
         g1p=0.d0 
         g2p=0.d0 
         g3p=0.d0 
         g4p=0.d0 
         g1a=0.d0 
         g2a=0.d0 
         g3a=0.d0 
         g4a=0.d0 
  166 pc=-pc-x2 
      g1m=g11-g1 
      g1mp=g1p1+g1p 
      gp=g2+g21 
      ge=g4+g41+x2*g1 
  167 go to(169,170),kw 
  155 do 171 k=1,5 
      z4=xxi(k)-1.d0 
      z1=exp(pc*z4) 
      y1=pc*xxi(k) 
      z2=1.d0+z1 
      z3=x2+y1 
      y2=pc*aai(k)*sqrt(pc*z3)/z2 
      y4=cci(k)+pc 
      z=y4+x2 
      y6=bbi(k)*sqrt(y4*z) 
      wo=wo+y2+y6 
      if((lst.eq.0).and.(kw.ne.1)) go to 172 
      z5=1./pc 
         y3=0.5d0*xxi(k)/z3-z4*z1/z2+1.5d0*z5 
      z6=1.d0/y4 
         y5=0.5d0*(1.d0/z+z6) 
      wop=wop+y2*y3+y6*y5 
      z3=y2/z3 
         z=y6/z 
  172 y2=y2*y1 
         y6=y6*y4 
      w1=w1+y2+y6 
      y3=y3+z5 
        y5=y5+z6 
      w1p=w1p+y2*y3+y6*y5 
      if(kw.eq.1) go to 171 
      if(lst.eq.0) go to 173 
      woa=woa+z3+z 
      z3=z3*y1 
         z=z*y4 
      w1a=w1a+z3+z 
  173 y2=y2*y1 
         y6=y6*y4 
      w2=w2+y2+y6 
      if(lst.eq.0) go to 171 
      w2p=w2p+y2*(y3+z5)+y6*(y5+z6) 
      w2a=w2a+z3*y1+z*y4 
  171 continue 
      go to 163 
  156 if(kk1.eq.0) go to 191 
      do 197 k=1,5 
      wk6(k)=hhi(k) 
         wk7(k)=ggi(k) 
      wk4(k)=sqrt(vvi(k)+x2) 
  197 wk5(k)=sqrt(zzi(k)+x2) 
      do 195 k=1,24 
  195 abcg(k)=0.d0 
      k1=0 
      do 198 i=1,3 
      x=i 
         x1=x-0.5d0 
      do 190 k=1,5 
      k2=k1+k 
        k3=k2+65 
      z=(x+ggsi(k))/pc2 
         y=x1/zzi(k) 
      z1=wk7(k)*(z-cpp(2))*cpp(4) 
         z2=wk6(k)*(y-cpp(2))*cpp(4) 
      z4=xxi(k)*wk7(k)/wk4(k) 
        z5=wk6(k)/wk5(k) 
         z6=cpp(5)*(z4+z5) 
      asp(i)=asp(i)+abc(k2)*uac(k3)+z1*wk4(k)+z2*wk5(k)+z6*cpp(1) 
      z8=cpp(5)*(z4/(vvi(k)+x2)+z5/(zzi(k)+x2)) 
      aspa(i)=aspa(i)+abc(k2)/uac(k3)+z1/wk4(k)+z2/wk5(k)-z8*cpp(1) 
      k4=k2+15 
      z1=wk7(k)*(cpp(3)-z)*cpp(1) 
         z2=wk6(k)*(cpp(3)-y)*cpp(1) 
      bsp(i)=bsp(i)+abc(k4)*uac(k3)+z1*wk4(k)+z2*wk5(k)-z6 
      bspa(i)=bspa(i)+abc(k4)/uac(k3)+z1/wk4(k)+z2/wk5(k)+z8 
      k4=k4+15 
      csp(i)=csp(i)+abc(k4)*uac(k3) 
      cspa(i)=cspa(i)+abc(k4)/uac(k3) 
      k4=k4+15 
      gsp(i)=gsp(i)+abc(k4)*uac(k3) 
      gspa(i)=gspa(i)+abc(k4)/uac(k3) 
      wk6(k)=wk6(k)*zzi(k) 
  190 wk7(k)=wk7(k)*vvi(k) 
  198 k1=k1+5 
      kk1=0 
  191 z=pc-pc1 
         z1=2.d0*z 
         z2=1.5d0*z 
      wo=gsp(1)+z*(csp(1)+z*(bsp(1)+z*asp(1))) 
      w1=gsp(2)+z*(csp(2)+z*(bsp(2)+z*asp(2))) 
      wop=csp(1)+z1*(bsp(1)+z2*asp(1)) 
      w1p=csp(2)+z1*(bsp(2)+z2*asp(2)) 
      if(kw.eq.1) go to 163 
      w2=gsp(3)+z*(csp(3)+z*(bsp(3)+z*asp(3))) 
      if(lst.eq.0) go to 163 
      w2p=csp(3)+z1*(bsp(3)+z2*asp(3)) 
      woa=gspa(1)+z*(cspa(1)+z*(bspa(1)+z*aspa(1))) 
      w1a=gspa(2)+z*(cspa(2)+z*(bspa(2)+z*aspa(2))) 
      w2a=gspa(3)+z*(cspa(3)+z*(bspa(3)+z*aspa(3))) 
      go to 163 
c 
c *** relativistic asymptotics (nz=4) 
    4 nzr=4 
  520  ro=pl*cu(58) 
          hi=1.d0/emue 
          r1=ro*0.5d0*hi 
          pit=pi2*al1 
          pt2=pit*al1 
          pa=pt2-1.5d0 
      hu=psi*al1+1.d0 
      do 525 it=1,4 
          hu1=hu 
          hu2=hu**2 
          hu =2.d0*(hu2*hu+r1)/(3.d0*hu2+pa) 
      if(abs(hu1/hu-1.d0).le.eit) go to 527 
  525  continue 
          r=r1**2+(pa*.3333333333d0)**3 
      x=(r1+sqrt(r))**(1.d0/3.d0) 
          hu=x-pa/(3.d0*x) 
  527  continue 
          hu2=hu**2 
          psi=-7.77d2 
          if(al1.gt.1.d-8) psi=(hu-1.d0)/al1 
          pe=.25d0*(hu2**2+2.d0*pa*hu2+.46666666667d0*pt2**2-pt2) 
          ee=(3.d0*pe+0.5d0*(3.d0*hu2+pt2))/ro-hi 
          pe=pe+1.1550858d0 
          ee=ee-1.1550858d0/ro 
          if(lst.eq.0) go to 555 
          r=1.d0/(3.d0*hu2+pa) 
          r1=hu2+pa 
          pt=pit*(hu2-0.5d0+0.46666666667d0*pt2)-2.d0*pit*hu2*r1*r 
          r2=hi*hu*r 
          ppl=r2*r1 
          et=pit*(hu2*(1.d0-4.d0*pt2*r)+1.4d0*pt2-0.5d0)/ro 
          epl=(3.d0*(ppl+r2)-ee-hi) 
  555  ee=ee*cu(59) 
       pe=pe*cu(60) 
      if(kpar.ne.1) go to 558 
      hpr=0.d0 
      if(t.lt.5.96d0) go to 558 
      eta=alf*hu 
      if(eta.gt.6.d1) go to 558 
      hu1=exp(-eta) 
      al2=al1**2 
      hpr=hu1*(1.2d1*al2-3.d0+hu1*((0.444444444d0*al2-1.d0) 
     * *hu1-1.5d0*(al2-1.d0)))/(eta*r1) 
  558 continue 
          if(lst.eq.0)go to 556 
          pt=pt*cu(61) 
          ppl=ppl*cu(59) 
          epl=epl*cu(59) 
          et=et*cu(2) 
 556   if(kentr.eq.0) go to 557 
          y=cu(62) 
      se=y*al1*(hu2+.466666666667d0*pt2-.5d0)/pl 
          if(lst.eq.0) go to 557 
      spl=-se+2.d0*pi2*cu(2)*al1*r2 
          st=se/t+2.d0*y*pt2*(0.46666666667d0-2.d0*hu2*r)/(pl*cu(1)) 
  557  go to 135 
c 
c *** interpolation between perfect gas and expansion over 
c *** half-integer f-d functions 
   52 pl1=x*emue 
         pl2=y*emue 
        nzp1=1 
         nzp2=2 
c 
c *** interpolation over density 
   83 kin=ki 
      lst1=lst 
         lst=1 
         kk=0 
      dni=pl 
      if(lst1.eq.0) go to 81 
      kk=1 
      tni=t 
         t=t*(cu(4)+dst) 
   81 pl=pl1 
         nz=nzp1 
        ki=2 
      go to 90 
   77 pn1=pe 
         en1=ee 
        sn1=se 
         psn1=psi 
         hprn1=hpr 
      pnp1=ppl 
         enp1=epl/pl 
         snp1=spl/pl 
      pl=pl2 
         nz=nzp2 
        ki=3 
         nzr1=nzr 
      go to 90 
   78 if(kk.eq.2) go to 92 
      wv=pl2-pl1 
         wv1=cu(4)/wv 
      wv4=cu(15)*wv1 
         wv3=dni-pl1 
         wv5=wv1*wv3 
   92 x=epl/pl 
      x1=pe-pn1 
         x2=x1*wv4 
         x3=x1*wv1 
      z1=pnp1+ppl-x2 
         z2=x3-z1-pnp1 
         z1=wv1*z1 
      pe=pn1+wv3*(pnp1+wv5*(z2+wv3*z1)) 
      x1=ee-en1 
         x2=x1*wv4 
         x3=x1*wv1 
      y1=enp1+x-x2 
         y2=x3-y1-enp1 
         y1=wv1*y1 
      ee=en1+wv3*(enp1+wv5*(y2+wv3*y1)) 
      if(kentr.eq.0) go to 91 
      x1=se-sn1 
         x2=x1*wv4 
         x3=x1*wv1 
      v1=snp1+spl/pl-x2 
         v2=x3-v1-snp1 
         v1=wv1*v1 
      se=sn1+wv3*(snp1+wv5*(v2+wv3*v1)) 
   91 if(kk.eq.1) go to 82 
      if(kk.eq.2) go to 124 
  127 x1=(dni-pl1)*wv1 
      psi=(psi-psn1)*x1+psn1 
      hpr=(hpr-hprn1)*x1+hprn1 
      nzr=10*nzr1+nzr 
       pl=dni 
        plm=pl/emue 
      ki=kin 
         lst=lst1 
  134 nz=0 
         pi=ei*pl 
  135  go to(57,77,78,79,80,577,578), ki 
   82 pn2=pe 
         en2=ee 
        sn2=se 
      kk=2 
         t=tni 
      go to 81 
  124 x1=t*dst 
      pt=(pn2-pe)/x1 
         et=(en2-ee)/x1 
      if(kentr.eq.0) go to 126 
      st=(sn2-se)/x1 
      spl=snp1+wv5*(cu(15)*v2+cu(17)*wv3*v1) 
      spl=dni*spl 
  126 ppl=pnp1+wv5*(cu(15)*z2+cu(17)*wv3*z1) 
      epl=enp1+wv5*(cu(15)*y2+cu(17)*wv3*y1) 
      epl=dni*epl 
      go to 127 
c 
c *** interpolation between degenerate gas and expansion over 
c *** half-integer fermi-dirac functions 
   62 pl1=x*emue 
         pl2=y*emue 
        nzp1=2 
         nzp2=3 
      go to 83 
c 
c *** interpolation over temperature 
  128 kit=ki 
      lst2=lst 
         lst=1 
         kkt=0 
      dnt=t 
      if(lst2.eq.0) go to 129 
      kkt=1 
         plni=pl 
        pl=pl*(cu(4)+dst) 
  129 t=t1 
         nz=nz1 
         ki=4 
      go to 90 
   79 pnt=pe 
         ent=ee 
         snt=se 
         psnt=psi 
         hprnt=hpr 
         pntt=pt 
         entt=et 
         sntt=st 
         t=t2 
         nz=nz2 
         ki=5 
         nzrt=nzr 
      go to 90 
   80 if(kkt.eq.2) go to 93 
      vw=t2-t1 
         vw1=cu(4)/vw 
      vw4=cu(15)*vw1 
         vw3=dnt-t1 
        vw5=vw1*vw3 
   93 x1=pe-pnt 
         x2=x1*vw4 
         x3=x1*vw1 
      z1=pntt+pt-x2 
         z2=x3-z1-pntt 
         z1=vw1*z1 
      pe=pnt+vw3*(pntt+vw5*(z2+vw3*z1)) 
      x1=ee-ent 
         x2=x1*vw4 
         x3=x1*vw1 
      y1=entt+et-x2 
         y2=x3-y1-entt 
         y1=vw1*y1 
      ee=ent+vw3*(entt+vw5*(y2+vw3*y1)) 
      if(kentr.eq.0) go to 94 
      x1=se-snt 
         x2=x1*vw4 
         x3=x1*vw1 
      v1=sntt+st-x2 
         v2=x3-v1-sntt 
         v1=vw1*v1 
      se=snt+vw3*(sntt+vw5*(v2+vw3*v1)) 
   94 if(kkt.eq.1) go to 130 
      if(kkt.eq.2) go to 131 
  133 x1=(dnt-t1)*vw1 
      psi=(psi-psnt)*x1+psnt 
      hpr=(hpr-hprnt)*x1+hprnt 
      nzr=10*nzrt+nzr 
      t=dnt 
         alf=cu(1)/t 
         al1=cu(4)/alf 
         ei=t*rg 
      eg=cu(3)*ei 
         ki=kit 
         lst=lst2 
      go to 134 
  130 pnt2=pe 
        ent2=ee 
         snt2=se 
      kkt=2 
         pl=plni 
      go to 129 
  131 x1=cu(4)/dst 
      ppl=x1*(pnt2-pe)/pl 
         epl=x1*(ent2-ee) 
      if(kentr.eq.0) go to 132 
      spl=(snt2-se)*x1 
      st=sntt+vw5*(cu(15)*v2+cu(17)*vw3*v1) 
  132 pt=pntt+vw5*(cu(15)*z2+cu(17)*vw3*z1) 
      et=entt+vw5*(cu(15)*y2+cu(17)*vw3*y1) 
      go to 133 


c..finish the pipeline loop
6969  continue
      enddo
      end 

