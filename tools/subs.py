subs = {
  'all' : [ ## Special key: Changes are applied to all applicable conversions automatically
    [None,None]
  ],
  'mixed' : [
    ['zc','ds'],
    ('ZC','DS'),
    ('zc','ds'),
    ('z_check','d_check'),
    ('complex\(kind=c_double_complex\)', 'real\(kind=c_double\)' ),
    ('real\(kind=c_float_complex\)',     'real\(kind=c_float\)'  ),
    ('PLASMA_Complex64_t','double'),
    ('PLASMA_Complex32_t','float'),
    ('PlasmaComplexDouble','PlasmaRealDouble'),
    ('PlasmaComplexFloat','PlasmaRealFloat'),
    ('zlange','dlange'),
    ('zlag2c','dlag2s'),
    ('clag2z','slag2d'),
    ('zlacpy','dlacpy'),
    ('zherfb','dsyrfb'),
    ('zherf','dsyrf'),
    ('zgemm','dgemm'),
    ('zherk','dsyrk'),
    ('zher2k','dsyr2k'),
    ('zlansy','dlansy'),
    ('zaxpy','daxpy'),
    ('zgeadd','dgeadd'),
    ('zplghe','dplgsy'),
    ('pzgetrf','pdgetrf'),
    ('pcgetrf','psgetrf'),
    ('pzhetrd','pdsytrd'),
    ('pzhbrdt','pdsbrdt'),
    ('ztrsm','dtrsm'),
    ('ctrsm','strsm'),
    ('CBLAS_SADDR',''),
    ('zlarnv','dlarnv'),
    ('zgesv','dgesv'),
    ('zhemm','dsymm'),
    ('zlanhe','dlansy'),
    ('zlaghe','dlagsy'),
    ('ztrmm','dtrmm'),
    ('ctrmm','strmm'),
    ('Conj',''),
    ('zpotrf','dpotrf'),
    ('cpotrf','spotrf'),
    ('zhegv','dsygv'),
    ('zhegst','dsygst'),
    ('PLASMA_Alloc_Workspace_zgels','PLASMA_Alloc_Workspace_dgels'),
    ('plasma_pc',    'plasma_ps'),
    ('plasma_pz',    'plasma_pd'),
    ('PLASMA_z',     'PLASMA_d'),
    ('PLASMA_c',     'PLASMA_s'),
    ('plasma_coop',  'plasma_soop'),
    ('plasma_zoop',  'plasma_doop'),
    ('plasma_cip',   'plasma_sip'),
    ('plasma_zip',   'plasma_dip'),
    ('plasma_cdesc', 'plasma_sdesc'),
    ('plasma_zdesc', 'plasma_ddesc'),
    ('unmqr','ormqr'),
    ('unmlq','ormlq'),
    ('he2hb','sy2sb'),
    ('hb2st','sy2st'),
  ],
  'normal' : [ ## Dictionary is keyed on substitution type
    ['s','d','c','z'], ## Special Line Indicating type columns

    ('#define PRECISION_s','#define PRECISION_d','#define PRECISION_c','#define PRECISION_z'),
    ('#undef PRECISION_s', '#undef PRECISION_d', '#undef PRECISION_c', '#undef PRECISION_z' ),
    ('#define REAL',  '#define REAL',  '#define COMPLEX','#define COMPLEX'),
    ('#undef COMPLEX','#undef COMPLEX','#undef REAL',    '#undef REAL'),
    ('#define SINGLE','#define DOUBLE','#define SINGLE', '#define DOUBLE'),
    ('#undef DOUBLE', '#undef SINGLE', '#undef DOUBLE',  '#undef SINGLE' ),

    # C and Fortran types
    ('real\(kind=c_float\)', 'real\(kind=c_double\)', 'complex\(kind=c_float_complex\)', 'complex\(kind=c_double_complex\)'),
    ('real\(kind=c_float\)', 'real\(kind=c_double\)', 'real\(kind=c_float\)', 'real\(kind=c_double\)'),
    ('real',           'double precision','real',                  'double precision'      ),
    ('REAL',           'DOUBLE_PRECISION','COMPLEX',               'COMPLEX_16'            ),
    ('float',          'double',          'float _Complex',        'double _Complex'       ),
    ('float',          'double',          'float',                 'double'                ),
    ('float',          'double',          'PLASMA_Complex32_t',    'PLASMA_Complex64_t'    ),
    ('float',          'double',          'PLASMA_voidComplex32_t','PLASMA_voidComplex64_t'),
    ('PlasmaRealFloat','PlasmaRealDouble','PlasmaComplexFloat',    'PlasmaComplexDouble'   ),
    ('real',           'double precision','complex',               'complex\(kind=wp\)'    ),
    ('sstedc','dstedc','cstedc','zstedc'),
    ('pssbcpy','pdsbcpy','pchbcpy','pzhbcpy'),
    ('psgbcpy','pdgbcpy','pcgbcpy','pzgbcpy'),
    ('ssbtype','dsbtype','chbtype','zhbtype'),
    ('cblas_snrm2','cblas_dnrm2','cblas_scnrm2','cblas_dznrm2'),
    ('cblas_sasum','cblas_dasum','cblas_scasum','cblas_dzasum'),
    ('CORE_sasum','CORE_dasum','CORE_scasum','CORE_dzasum'),
    ('core_sasum','core_dasum','core_scasum','core_dzasum'),
    ('qwrapper_sasum','qwrapper_dasum','qwrapper_scasum','qwrapper_dzasum'),
    ('ipt_s','ipt_d','ipt_c','ipt_z'),
    ('coreblas_s','coreblas_d','coreblas_c','coreblas_z'),
    ('sytra1','sytra1','hetra1','hetra1'),
    ('sgecon','dgecon','cgecon','zgecon'),
    ('spocon','dpocon','cpocon','zpocon'),
    ('slacn2','dlacn2','clacn2','zlacn2'),
    ('SLACN2','DLACN2','CLACN2','ZLACN2'),
    ('ssygst','dsygst','chegst','zhegst'),
    ('SSYGST','DSYGST','CHEGST','ZHEGST'),
    ('ssterf','dsterf','ssterf','dsterf'),
    ('ssytrd','dsytrd','chetrd','zhetrd'),
    ('SSYTRD','DSYTRD','CHETRD','ZHETRD'),
    ('STILE','DTILE','CTILE','ZTILE'),
    ('stile','dtile','ctile','ztile'),
    ('slag2d','dlag2s','clag2z','zlag2c'),
    ('ssyrfb','dsyrfb','cherfb','zherfb'),
    ('ssyrf','dsyrf','cherf','zherf'),
    ('sbarrier','dbarrier','cbarrier','zbarrier'),
    ('saxpy','daxpy','caxpy','zaxpy'),
    ('sgeadd','dgeadd','cgeadd','zgeadd'),
    ('sgeam','dgeam','cgeam','zgeam'),
    ('ssymm','dsymm','chemm','zhemm'),
    ('SSYMM','DSYMM','CHEMM','ZHEMM'),
    ('ssymv','dsymv','chemv','zhemv'),
    ('SSYMV','DSYMV','CHEMV','ZHEMV'),
    ('ssyrk','dsyrk','cherk','zherk'),
    ('SSYRK','DSYRK','CHERK','ZHERK'),
    ('ssyr2','dsyr2','cher2','zher2'),
    ('SSYR2','DSYR2','CHER2','ZHER2'),
    ('sgesv','dgesv','cgesv','zgesv'),
    ('SUNGESV','SUNGESV','CUNGESV','CUNGESV'),
    ('SGESV','SGESV','CGESV','CGESV'),
    ('SGESV','DGESV','CGESV','ZGESV'),
    ('sgels','dgels','cgels','zgels'),
    ('SGELS','DGELS','CGELS','ZGELS'),
    ('ssyev','dsyev','cheev','zheev'),
    ('ssyevd','dsyevd','cheevd','zheevd'),
    ('ssyevr','dsyevr','cheevr','zheevr'),
    ('SSYEV','DSYEV','CHEEV','ZHEEV'),
    ('SSYEVD','DSYEVD','CHEEVD','ZHEEVD'),
    ('SSYEVR','DSYEVR','CHEEVR','ZHEEVR'),
    ('ssygv','dsygv','chegv','zhegv'),
    ('SSYGV','DSYGV','CHEGV','ZHEGV'),
    ('sgemm','dgemm','cgemm','zgemm'),
    ('SGEMM','DGEMM','CGEMM','ZGEMM'),
    ('sposv','dposv','cposv','zposv'),
    ('SPOSV','SPOSV','CPOSV','CPOSV'),
    ('SPOSV','DPOSV','CPOSV','ZPOSV'),
    ('SGEBRD','DGEBRD','CGEBRD','ZGEBRD'),
    ('sgebrd','dgebrd','cgebrd','zgebrd'),
    ('sgeev','dgeev','cgeev','zgeev'),
    ('SGEEV','DGEEV','CGEEV','ZGEEV'),
    ('SGEHRD','DGEHRD','CGEHRD','ZGEHRD'),
    ('sgehrd','dgehrd','cgehrd','zgehrd'),
    ('sgerbb','dgerbb','cgerbb','zgerbb'),
    ('sgerbh','dgerbh','cgerbh','zgerbh'),
    ('sgerbbrh','dgerbbrh','cgerbbrh','zgerbbrh'),
    ('sgbrdb','dgbrdb','cgbrdb','zgbrdb'),
    ('SGESVD','DGESVD','CGESVD','ZGESVD'),
    ('SGESDD','DGESDD','CGESDD','ZGESDD'),
    ('sgesvd','dgesvd','cgesvd','zgesvd'),
    ('sgesdd','dgesdd','cgesdd','zgesdd'),
    ('sbdsdc','dbdsdc','sbdsdc','dbdsdc'),
    ('ssymm','dsymm','csymm','zsymm'),
    ('SSYMM','DSYMM','CSYMM','ZSYMM'),
    ('ssyrk','dsyrk','csyrk','zsyrk'),
    ('SSYRK','DSYRK','CSYRK','ZSYRK'),
    ('ssyr2k','dsyr2k','csyr2k','zsyr2k'),
    ('SSYR2K','DSYR2K','CSYR2K','ZSYR2K'),
    ('strgmm','dtrgmm','ctrgmm','ztrgmm'),
    ('strbmm','dtrbmm','ctrbmm','ztrbmm'),
    ('strmm','dtrmm','ctrmm','ztrmm'),
    ('STRMM','DTRMM','CTRMM','ZTRMM'),
    ('strsm','dtrsm','ctrsm','ztrsm'),
    ('STRSM','DTRSM','CTRSM','ZTRSM'),
    ('sgelq3','dgelq3','cgelq3','zgelq3'),
    ('sgelq2','dgelq2','cgelq2','zgelq2'),
    ('sgelqf','dgelqf','cgelqf','zgelqf'),
    ('sgelqfrh','dgelqfrh','cgelqfrh','zgelqfrh'),
    ('SGELQF','DGELQF','CGELQF','ZGELQF'),
    ('sgelqs','dgelqs','cgelqs','zgelqs'),
    ('SGELQS','DGELQS','CGELQS','ZGELQS'),
    ('sgeqr2','dgeqr2','cgeqr2','zgeqr2'),
    ('sgeqr3','dgeqr3','cgeqr3','zgeqr3'),
    ('sgeqp3','dgeqp3','cgeqp3','zgeqp3'),
    ('sgeqrf','dgeqrf','cgeqrf','zgeqrf'),
    ('sgeqrfrh','dgeqrfrh','cgeqrfrh','zgeqrfrh'),
    ('SGEQRF','DGEQRF','CGEQRF','ZGEQRF'),
    ('sgeqrs','dgeqrs','cgeqrs','zgeqrs'),
    ('SGEQRS','DGEQRS','CGEQRS','ZGEQRS'),
    ('sgetf2','dgetf2','cgetf2','zgetf2'),
    ('sgetrf','dgetrf','cgetrf','zgetrf'),
    ('SGETRF','DGETRF','CGETRF','ZGETRF'),
    ('sgetrs','dgetrs','cgetrs','zgetrs'),
    ('SGETRS','DGETRS','CGETRS','ZGETRS'),
    ('slacgv','dlacgv','clacgv','zlacgv'),
    ('slacpy','dlacpy','clacpy','zlacpy'),
    ('SLACPY','DLACPY','CLACPY','ZLACPY'),
    ('slagsy','dlagsy','claghe','zlaghe'),
    ('slagsy','dlagsy','clagsy','zlagsy'),
    ('SLANGE','DLANGE','CLANGE','ZLANGE'),
    ('SLANSY','DLANSY','CLANHE','ZLANHE'),
    ('SLANSY','DLANSY','CLANSY','ZLANSY'),
    ('SLANTR','DLANTR','CLANTR','ZLANTR'),
    ('slaset','dlaset','claset','zlaset'),
    ('SLASET','DLASET','CLASET','ZLASET'),
    ('slange','dlange','clange','zlange'),
    ('slansy','dlansy','clanhe','zlanhe'),
    ('slansy','dlansy','clansy','zlansy'),
    ('slantr','dlantr','clantr','zlantr'),
    ('slarfb','dlarfb','clarfb','zlarfb'),
    ('slarfg','dlarfg','clarfg','zlarfg'),
    ('slarft','dlarft','clarft','zlarft'),
    ('slarfx','dlarfx','clarfx','zlarfx'),
    ('slarfy','dlarfy','clarfy','zlarfy'),
    ('SLARFY','DLARFY','CLARFY','ZLARFY'),
    ('slarfw','dlarfw','clarfw','zlarfw'),
    ('slarnv','dlarnv','clarnv','zlarnv'),
    ('slaswp','dlaswp','claswp','zlaswp'),
    ('SLASWP','DLASWP','CLASWP','ZLASWP'),
    ('slatms','dlatms','clatms','zlatms'),
    ('slatro','dlatro','clatro','zlatro'),
    ('SPEMV','DPEMV','CPEMV','ZPEMV'),
    ('splgsy','dplgsy','cplghe','zplghe'),
    ('spotrf','dpotrf','cpotrf','zpotrf'),
    ('spotrf','dpotrf','cpotrf','zpotrf'),
    ('SPOTRF','DPOTRF','CPOTRF','ZPOTRF'),
    ('spotrs','dpotrs','cpotrs','zpotrs'),
    ('SPOTRS','DPOTRS','CPOTRS','ZPOTRS'),
    ('sorgbr','dorgbr','cungbr','zungbr'),
    ('sorgbrrh','dorgbrrh','cungbrrh','zungbrrh'),
    ('SORGBR','DORGBR','CUNGBR','ZUNGBR'),
    ('sorghr','dorghr','cunghr','zunghr'),
    ('SORGHR','DORGHR','CUNGHR','ZUNGHR'),
    ('sorglq','dorglq','cunglq','zunglq'),
    ('sorglqrh','dorglqrh','cunglqrh','zunglqrh'),
    ('SORGLQ','DORGLQ','CUNGLQ','ZUNGLQ'),
    ('sorgqr','dorgqr','cungqr','zungqr'),
    ('sorgqrrh','dorgqrrh','cungqrrh','zungqrrh'),
    ('SORGQR','DORGQR','CUNGQR','ZUNGQR'),
    ('sorgtr','dorgtr','cungtr','zungtr'),
    ('SORGTR','DORGTR','CUNGTR','ZUNGTR'),
    ('sormlq','dormlq','cunmlq','zunmlq'),
    ('sormlqrh','dormlqrh','cunmlqrh','zunmlqrh'),
    ('SORMLQ','DORMLQ','CUNMLQ','ZUNMLQ'),
    ('sormqr','dormqr','cunmqr','zunmqr'),
    ('sormqrrh','dormqrrh','cunmqrrh','zunmqrrh'),
    ('SORMQR','DORMQR','CUNMQR','ZUNMQR'),
    ('ssytrd','dsytrd','chetrd','zhetrd'),
    ('ssbrdb','dsbrdb','chbrdb','zhbrdb'),
    ('SSBRDB','DSBRDB','CHBRDB','ZHBRDB'),
    ('ssbrdt','dsbrdt','chbrdt','zhbrdt'),
    ('SSBRDT','DSBRDT','CHBRDT','ZHBRDT'),
    ('slamch','dlamch','slamch','dlamch'),
    ('slarnv','dlarnv','slarnv','dlarnv'),
    ('slauum','dlauum','clauum','zlauum'),
    ('SLAUUM','DLAUUM','CLAUUM','ZLAUUM'),
    ('spotri','dpotri','cpotri','zpotri'),
    ('SPOTRI','DPOTRI','CPOTRI','ZPOTRI'),
    ('strtri','dtrtri','ctrtri','ztrtri'),
    ('STRTRI','DTRTRI','CTRTRI','ZTRTRI'),
    ('sgetri','dgetri','cgetri','zgetri'),
    ('SGETRI','DGETRI','CGETRI','ZGETRI'),
    ('sshift','dshift','cshift','zshift'),
    ('sgetmo','dgetmo','cgetmo','zgetmo'),
    ('sgetmi','dgetmi','cgetmi','zgetmi'),
    ('SGETMI','DGETMI','CGETMI','ZGETMI'),
    ('sgecfi','dgecfi','cgecfi','zgecfi'),
    ('SGECFI','DGECFI','CGECFI','ZGECFI'),
    ('spack','dpack','cpack','zpack'),
    ('strsmpl','dtrsmpl','ctrsmpl','ztrsmpl'),
    ('STRSMPL','DTRSMPL','CTRSMPL','ZTRSMPL'),
    ('splgsy','dplgsy','cplghe','zplghe'),
    ('SPLGSY','DPLGSY','CPLGHE','ZPLGHE'),
    ('splgsy','dplgsy','cplgsy','zplgsy'),
    ('SPLGSY','DPLGSY','CPLGSY','ZPLGSY'),
    ('splrnt','dplrnt','cplrnt','zplrnt'),
    ('SPLRNT','DPLRNT','CPLRNT','ZPLRNT'),
    ('spltmg','dpltmg','cpltmg','zpltmg'),
    ('SPLTMG','DPLTMG','CPLTMG','ZPLTMG'),
    ('\*\*T','\*\*T','\*\*H','\*\*H'),
    ('BLAS_s','BLAS_d','BLAS_s','BLAS_d'),
    ('BLAS_s','BLAS_d','BLAS_c','BLAS_z'),
    ('fabsf','fabs','cabsf','cabs'),
    ('imagf','imag','imagf','imag'),
    ('cosf','cos','ccosf','ccos'),
    ('powf','pow','cpowf','cpow'),
    ('cblas_sscal','cblas_dscal','cblas_csscal','cblas_zdscal'),
    ('cblas_is','cblas_id','cblas_ic','cblas_iz'),
    ('cblas_is','cblas_id','cblas_is','cblas_id'),
    ('cblas_s','cblas_d','cblas_c','cblas_z'),
    ('','','CBLAS_SADDR','CBLAS_SADDR'),
    ('CblasTrans','CblasTrans','CblasConjTrans','CblasConjTrans'),
    ('','','conjf','conj'),
    ('CORE_S','CORE_D','CORE_C','CORE_Z'),
    ('CORE_s','CORE_d','CORE_c','CORE_z'),
    ('CORE_s','CORE_d','CORE_s','CORE_d'),
    ('core_s','core_d','core_c','core_z'),
    ('qwrapper_s','qwrapper_d','qwrapper_c','qwrapper_z'),
    ('example_s','example_d','example_c','example_z'),
    ('lapack_s','lapack_d','lapack_c','lapack_z'),
    ('lapack_slamch','lapack_dlamch','lapack_slamch','lapack_dlamch'),
    ('plasma_ps','plasma_pd','plasma_pc','plasma_pz'),
    ('PLASMA_s','PLASMA_d','PLASMA_c','PLASMA_z'),
    ('PLASMA_S','PLASMA_D','PLASMA_C','PLASMA_Z'),
    ('plasma_s','plasma_d','plasma_c','plasma_z'),
    ('control_s','control_d','control_c','control_z'),
    ('compute_s','compute_d','compute_c','compute_z'),
    ('PLASMA_sor','PLASMA_dor','PLASMA_cun','PLASMA_zun'),
    ('PlasmaTrans','PlasmaTrans','PlasmaConjTrans','PlasmaConjTrans'),
    ('testing_ds','testing_ds','testing_zc','testing_zc'),
    ('time_s','time_d','time_c','time_z'),
    ('testing_s','testing_d','testing_c','testing_z'),
    ('TESTING_S','TESTING_D','TESTING_C','TESTING_Z'),
    ('stesting','dtesting','ctesting','ztesting'),
    ('SAUXILIARY','DAUXILIARY','CAUXILIARY','ZAUXILIARY'),
    ('sauxiliary','dauxiliary','cauxiliary','zauxiliary'),
    ('s_check','d_check','c_check','z_check'),
    ('ger','ger','gerc','gerc'),
    ('ger','ger','geru','geru'),
    ('symm','symm','hemm','hemm'),
    ('syrk','syrk','herk','herk'),
    ('syr2k','syr2k','her2k','her2k'),
    ('lansy','lansy','lanhe','lanhe'),
    ('plgsy','plgsy','plghe','plghe'),
    ('ormqr','ormqr','unmqr','unmqr'),
    ('','','crealf','creal'),
    ('slasrt_','dlasrt_','slasrt_','dlasrt_'),
    ('LAPACKE_s','LAPACKE_d','LAPACKE_c','LAPACKE_z'),
    ('SLAPACK','DLAPACK','CLAPACK','ZLAPACK'),
    ('slapack','dlapack','clapack','zlapack'),
    ('workspace_s','workspace_d','workspace_c','workspace_z'),
    ('Workspace_s','Workspace_d','Workspace_c','Workspace_z'),
    ('slassq', 'dlassq', 'classq', 'zlassq'),
    ('SLASSQ', 'DLASSQ', 'CLASSQ', 'ZLASSQ'),
  ],
  'tracing' : [
    ['plain','tau'],
    ('(\w+\*?)\s+(\w+)\s*\(([a-z* ,A-Z_0-9]*)\)\s*{\s+(.*)\s*#pragma tracing_start\s+(.*)\s+#pragma tracing_end\s+(.*)\s+}',r'\1 \2(\3){\n\4tau("\2");\5tau();\6}'),
    ('\.c','.c.tau'),
  ],
};
