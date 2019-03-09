use blasrs::level2::complex;
use blasrs::unstable::matrix;
use num_complex::Complex;
mod fixtures;

#[macro_use]
mod utils;

#[test]
fn gbmv() {
    let mut y = fixtures::complex::vector(6);
    complex::gbmv(
        'n',
        6,
        8,
        3,
        2,
        Complex::new(0.0, 0.0),
        &fixtures::complex::bandmatrix_nxm_ku_kl(8, 6, 3, 2),
        6,
        &fixtures::complex::vector(8),
        1,
        Complex::new(2.5, 0.5),
        &mut y,
        1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(-2.0085962414741516, 1.6058503985404968),
            Complex::new(-0.18816410377621651, -0.60837343707680702),
            Complex::new(1.8727443814277649, -0.72995787858963013),
            Complex::new(2.9619127362966537, -0.49696569144725800),
            Complex::new(-0.13906472176313400, 2.5643529072403908),
            Complex::new(-0.15649497509002686, -0.74832186102867126)
        ]
    );

    let mut y = fixtures::complex::vector(6);
    complex::gbmv(
        'n',
        6,
        8,
        3,
        2,
        Complex::new(0.2, 0.8),
        &fixtures::complex::bandmatrix_nxm_ku_kl(8, 6, 3, 2),
        6,
        &fixtures::complex::vector(8),
        1,
        Complex::new(2.5, 0.5),
        &mut y,
        1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(-1.0943454405994295, 1.0642741157616680),
            Complex::new(-0.38628598141225878, -1.2014790478399264),
            Complex::new(0.88995188622798471, -1.1623543999107093),
            Complex::new(7.9118205360687366E-002, 0.76470894087791352),
            Complex::new(-5.7862862678327798E-002, 0.62930340373117366),
            Complex::new(-0.36689070841052379, -0.18794836059332087)
        ]
    );

    let mut y = fixtures::complex::vector(8);
    complex::gbmv(
        't',
        6,
        8,
        3,
        2,
        Complex::new(0.2, 0.8),
        &fixtures::complex::bandmatrix_nxm_ku_kl(8, 6, 3, 2),
        6,
        &fixtures::complex::vector(6),
        1,
        Complex::new(0.0, 0.0),
        &mut y,
        1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(-0.38810574366012285, -1.6303050664096645),
            Complex::new(2.3103061857743366, -0.97417280945264095),
            Complex::new(-4.2076520090095242, -0.50591381502069677),
            Complex::new(-0.81853765257644051, 1.4503261671516843),
            Complex::new(2.6631576242668459E-002, -0.63622673851522416),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
        ]
    );

    let mut y = fixtures::complex::vector(8);
    complex::gbmv(
        'c',
        6,
        8,
        3,
        2,
        Complex::new(0.2, 0.8),
        &fixtures::complex::bandmatrix_nxm_ku_kl(8, 6, 3, 2),
        6,
        &fixtures::complex::vector(6),
        -1,
        Complex::new(1.0, 0.0),
        &mut y,
        -1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(-0.64901006221771240, 0.77214217185974121),
            Complex::new(-0.11916876584291458, -0.21951562166213989),
            Complex::new(0.66413569450378418, -0.42481029033660889),
            Complex::new(0.76170726548233147, -0.74829467121041837),
            Complex::new(0.49362246616457550, -0.11820536025789097),
            Complex::new(-0.14226157594261568, -3.1888940811820188),
            Complex::new(-2.0779711022358427, 1.5114922978334300),
            Complex::new(-2.4265906255376652, 0.38932194147661164),
        ]
    );

    let mut y = fixtures::complex::vector(8);
    complex::gbmv(
        'c',
        6,
        8,
        3,
        2,
        Complex::new(0.0, 0.0),
        &fixtures::complex::bandmatrix_nxm_ku_kl(8, 6, 3, 2),
        6,
        &fixtures::complex::vector(6),
        -1,
        Complex::new(1.0, 0.0),
        &mut y,
        -1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(-0.64901006221771240, 0.77214217185974121),
            Complex::new(-0.11916876584291458, -0.21951562166213989),
            Complex::new(0.66413569450378418, -0.42481029033660889),
            Complex::new(1.1009690761566162, -0.41898009181022644),
            Complex::new(0.14377148449420929, 0.99698686599731445),
            Complex::new(-0.11775359511375427, -0.27577802538871765),
            Complex::new(-0.91206836700439453, 1.2560187578201294),
            Complex::new(-1.4375861883163452, 0.64667439460754395),
        ]
    );

    let result = std::panic::catch_unwind(|| {
        complex::gbmv(
            'x',
            6,
            8,
            4,
            4,
            Complex::new(0.0, 0.0),
            &vec![],
            6,
            &vec![],
            1,
            Complex::new(2.0, 0.5),
            &mut vec![],
            1,
        )
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::gbmv(
            't',
            -6,
            8,
            4,
            4,
            Complex::new(0.0, 0.0),
            &vec![],
            6,
            &vec![],
            1,
            Complex::new(2.0, 0.5),
            &mut vec![],
            1,
        )
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::gbmv(
            't',
            6,
            -8,
            4,
            4,
            Complex::new(0.0, 0.0),
            &vec![],
            6,
            &vec![],
            1,
            Complex::new(2.0, 0.5),
            &mut vec![],
            1,
        )
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::gbmv(
            't',
            6,
            8,
            -4,
            4,
            Complex::new(0.0, 0.0),
            &vec![],
            6,
            &vec![],
            1,
            Complex::new(2.0, 0.5),
            &mut vec![],
            1,
        )
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::gbmv(
            't',
            6,
            8,
            4,
            -4,
            Complex::new(0.0, 0.0),
            &vec![],
            6,
            &vec![],
            1,
            Complex::new(2.0, 0.5),
            &mut vec![],
            1,
        )
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::gbmv(
            't',
            6,
            8,
            4,
            4,
            Complex::new(0.0, 0.0),
            &vec![],
            6,
            &vec![],
            1,
            Complex::new(2.0, 0.5),
            &mut vec![],
            1,
        )
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::gbmv(
            't',
            6,
            8,
            4,
            4,
            Complex::new(0.0, 0.0),
            &vec![],
            6,
            &vec![],
            0,
            Complex::new(2.0, 0.5),
            &mut vec![],
            1,
        )
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::gbmv(
            't',
            6,
            8,
            4,
            4,
            Complex::new(0.0, 0.0),
            &vec![],
            6,
            &vec![],
            1,
            Complex::new(2.0, 0.5),
            &mut vec![],
            0,
        )
    });
    assert!(result.is_err());
}

#[test]
fn gemv() {
    let mut y = fixtures::complex::vector(6);
    complex::gemv(
        'n',
        6,
        8,
        Complex::new(0.8, 0.2),
        &fixtures::complex::matrix_mxn(6, 8),
        6,
        &fixtures::complex::vector(8),
        1,
        Complex::new(0.2, 0.8),
        &mut y,
        1,
    );
    println!("{:?}", y);
    capproximately!(
        y,
        vec![
            Complex::new(-4.0955083795571205, 1.1357517892058278),
            Complex::new(0.93896042218602016, 1.0838676654345458),
            Complex::new(-1.5853445176637875, 1.6962449347899724),
            Complex::new(-0.25555151836345036, -4.3800522220206233),
            Complex::new(-1.2051793139138998, -2.9157429210651786),
            Complex::new(0.38280021798680286, 0.19114175133433564),
        ]
    );

    let mut y = fixtures::complex::vector(6);
    complex::gemv(
        'n',
        6,
        8,
        Complex::new(0.8, 0.2),
        &fixtures::complex::matrix_mxn(6, 8),
        6,
        &fixtures::complex::vector(8),
        -1,
        Complex::new(0.0, 0.0),
        &mut y,
        -1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(2.6017477475494872, -0.62451989822044895),
            Complex::new(1.8074832963400977, -2.2773562366746307),
            Complex::new(-4.2649596122103990, 2.1071421255650922),
            Complex::new(0.27753542357045247, -1.9840473842815003),
            Complex::new(2.2299686271484305, 0.24845760476943476),
            Complex::new(0.20069066517622025, -1.0213684617581498),
        ]
    );

    let mut y = fixtures::complex::vector(6);
    complex::gemv(
        'n',
        6,
        8,
        Complex::new(0.8, 0.2),
        &fixtures::complex::matrix_mxn(6, 8),
        6,
        &fixtures::complex::vector(8),
        1,
        Complex::new(1.0, 0.0),
        &mut y,
        1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(-3.9970026807046444, 2.2726735819034305),
            Complex::new(0.66801290992029672, 1.0035901828539755),
            Complex::new(-1.3938842013734591, 0.82508814026657729),
            Complex::new(0.29003966083788602, -5.5960115682700131),
            Complex::new(-0.29257262206414292, -2.2331706205478370),
            Complex::new(6.7974918648228133E-002, 6.4722209339979742E-002),
        ]
    );

    let mut y = fixtures::complex::vector(8);
    complex::gemv(
        't',
        6,
        8,
        Complex::new(0.8, 0.2),
        &fixtures::complex::matrix_mxn(6, 8),
        6,
        &fixtures::complex::vector(6),
        1,
        Complex::new(1.0, 0.0),
        &mut y,
        1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(-1.0722266311585562, 1.4276441972981582),
            Complex::new(4.0697286775452355, -1.3751872740072097),
            Complex::new(-1.7667910641302869, -1.3146598668049647),
            Complex::new(2.7682588305027580, 0.43286241964465866),
            Complex::new(0.61287748489518667, -1.2227666045246042),
            Complex::new(-5.6604442768105556E-002, 1.3467393639005534),
            Complex::new(-2.2131934192050711, 4.1552655661274187),
            Complex::new(-2.9062411441513247, 2.3044710975088387),
        ]
    );

    let mut y = fixtures::complex::vector(8);
    complex::gemv(
        'c',
        6,
        8,
        Complex::new(0.8, 0.2),
        &fixtures::complex::matrix_mxn(6, 8),
        6,
        &fixtures::complex::vector(6),
        1,
        Complex::new(1.0, 0.0),
        &mut y,
        1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(2.2502367460895774, 2.0732672240604177),
            Complex::new(1.5386279370583802, 1.2491393341728396),
            Complex::new(3.5354757160882668, 0.53153277666431009),
            Complex::new(-0.65968123771844311, -0.51832467085111289),
            Complex::new(-8.1751634328745854E-002, 0.52362043122181512),
            Complex::new(-1.8050991363430615, -0.94441463917334822),
            Complex::new(-0.25656634156597757, 1.6792353267609732),
            Complex::new(-2.5932578406614151, -3.5422230487806061),
        ]
    );

    let mut y = fixtures::complex::vector(6);
    complex::gemv(
        'n',
        6,
        8,
        Complex::new(0.0, 0.0),
        &fixtures::complex::matrix_mxn(6, 8),
        6,
        &fixtures::complex::vector(8),
        1,
        Complex::new(0.2, 0.8),
        &mut y,
        1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(-0.74751576107018813, -0.36477962083786153),
            Complex::new(0.15177874642280853, -0.13923813908156957),
            Complex::new(0.47267537821345584, 0.44634650418678667),
            Complex::new(0.55537789695527984, 0.79697925443916384),
            Complex::new(-0.76883520735554689, 0.31441456547997237),
            Complex::new(0.19707170422482045, -0.14935848339436175),
        ]
    );

    let mut y = fixtures::complex::vector(6);
    let expect = fixtures::complex::vector(6);
    complex::gemv(
        'n',
        6,
        8,
        Complex::new(0.0, 0.0),
        &fixtures::complex::matrix_mxn(6, 8),
        6,
        &fixtures::complex::vector(6),
        1,
        Complex::new(1.0, 0.0),
        &mut y,
        1,
    );
    capproximately!(y, expect);

    let result = std::panic::catch_unwind(|| {
        complex::gemv(
            'x',
            6,
            8,
            Complex::new(0.0, 0.0),
            &vec![],
            6,
            &vec![],
            1,
            Complex::new(2.0, 0.5),
            &mut vec![],
            1,
        )
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::gemv(
            't',
            -6,
            8,
            Complex::new(0.0, 0.0),
            &vec![],
            6,
            &vec![],
            1,
            Complex::new(2.0, 0.5),
            &mut vec![],
            1,
        )
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::gemv(
            't',
            6,
            -8,
            Complex::new(0.0, 0.0),
            &vec![],
            6,
            &vec![],
            1,
            Complex::new(2.0, 0.5),
            &mut vec![],
            1,
        )
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::gemv(
            't',
            6,
            8,
            Complex::new(0.0, 0.0),
            &vec![],
            6,
            &vec![],
            1,
            Complex::new(2.0, 0.5),
            &mut vec![],
            1,
        )
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::gemv(
            't',
            6,
            8,
            Complex::new(0.0, 0.0),
            &vec![],
            6,
            &vec![],
            1,
            Complex::new(2.0, 0.5),
            &mut vec![],
            1,
        )
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::gemv(
            't',
            6,
            8,
            Complex::new(0.0, 0.0),
            &vec![],
            6,
            &vec![],
            1,
            Complex::new(2.0, 0.5),
            &mut vec![],
            1,
        )
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::gemv(
            't',
            6,
            8,
            Complex::new(0.0, 0.0),
            &vec![],
            6,
            &vec![],
            0,
            Complex::new(2.0, 0.5),
            &mut vec![],
            1,
        )
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::gemv(
            't',
            6,
            8,
            Complex::new(0.0, 0.0),
            &vec![],
            6,
            &vec![],
            1,
            Complex::new(2.0, 0.5),
            &mut vec![],
            0,
        )
    });
    assert!(result.is_err());
}

#[test]
fn gerc() {
    let mut a = fixtures::complex::matrix_mxn(6, 8);
    let y = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
        Complex::new(-0.9120683669483379, 1.2560188173061),
        Complex::new(-1.4375862408299789, 0.6466743904953449),
    ];

    complex::gerc(
        6,
        8,
        Complex::new(0.2, 0.8),
        &fixtures::complex::vector(6),
        1,
        &y,
        1,
        &mut a,
        6,
    );
    capproximately!(
        a,
        vec![
            Complex::new(1.4664377569938685, 1.8060944675544837),
            Complex::new(-0.53225092996643109, -0.45634091440603552),
            Complex::new(1.3676711769934335, 0.58364817261230650),
            Complex::new(1.5273647960307939, -1.2254245383182520),
            Complex::new(1.1563959711297374, 2.1474949691654932),
            Complex::new(-1.7831775159348064, 0.50551385861110387),
            Complex::new(-0.75941169595474800, -0.57340502179859332),
            Complex::new(-0.28224278033580208, -0.78213264706952679),
            Complex::new(-0.16007534447744504, -1.1160014767921838),
            Complex::new(2.1635202157291209, -1.0386515298872647),
            Complex::new(0.78619566928364970, -1.7700218301387962),
            Complex::new(-0.78970753495767643, 1.2175962397684426),
            Complex::new(-1.4891467992751752, 0.27223155051811232),
            Complex::new(-0.12951010544349723, -0.25532453277181810),
            Complex::new(-0.17490710560589115, 0.76336937172774166),
            Complex::new(-0.38122952832127499, 0.38852989041687136),
            Complex::new(-0.39195398537726411, 2.3235693548980874),
            Complex::new(-0.69758972915325101, -0.81081531972369780),
            Complex::new(0.43568331003189087, -5.4877474904060364E-002),
            Complex::new(-1.2375384569168091, 0.25014132261276245),
            Complex::new(-0.22426788508892059, 0.61824327707290649),
            Complex::new(0.37739565968513489, -0.17262350022792816),
            Complex::new(0.13333636522293091, -2.2239003181457520),
            Complex::new(0.80418950319290161, -1.2636144161224365),
            Complex::new(-0.52825871591104467, 1.0515473739625048),
            Complex::new(0.38661084813518909, -0.18238536895940438),
            Complex::new(1.5987282580693853, -1.3477283962554565),
            Complex::new(0.18347150310255600, -0.55494690365844279),
            Complex::new(-1.0816686910303539, -3.2462525018476240E-003),
            Complex::new(-7.3848983873964338E-002, 2.4312089369662865E-002),
            Complex::new(-4.7085681056201645E-002, -1.5882927279464436),
            Complex::new(-0.52236195772798288, 0.42419417151368499),
            Complex::new(-0.61206211258427001, 0.32620723072294699),
            Complex::new(-0.93465874966505846, 0.12460202594829736),
            Complex::new(0.73057521316282659, -0.22989490997416154),
            Complex::new(1.1698956214925516, 0.32927391868466860),
            Complex::new(1.2157758131530922, 2.5345520058413213),
            Complex::new(-0.74283140469563169, -0.38987560759970030),
            Complex::new(1.3678114595993076, 0.32901162588415334),
            Complex::new(0.21513199437223562, -0.15203527634896519),
            Complex::new(2.8540439631030279, 1.0935453027260844),
            Complex::new(0.19340614948484969, -1.6512506222549119),
            Complex::new(0.38594072116787248, 7.9234395057932172E-002),
            Complex::new(-1.1402800589724982, -0.19270504442970846),
            Complex::new(-1.5574412841395437, -0.95309580634241586),
            Complex::new(-1.3486081371207643, 0.89977827692739010),
            Complex::new(-0.25519137193589581, 0.80878144081759595),
            Complex::new(0.77664318903912266, -0.71173479570372455),
        ]
    );

    let mut a = fixtures::complex::matrix_mxn(6, 8);
    complex::gerc(
        6,
        8,
        Complex::new(0.2, 0.8),
        &fixtures::complex::vector(6),
        -1,
        &fixtures::complex::vector(8),
        -1,
        &mut a,
        6,
    );
    capproximately!(
        a,
        vec![
            Complex::new(0.88306036814678868, 1.0794348476983751),
            Complex::new(0.98235736678213703, -0.38432509127926684),
            Complex::new(1.0467817773918699, -0.26657089566477299),
            Complex::new(0.88155860714158418, -1.2266749209248053),
            Complex::new(0.10640467108044993, 1.8599184957333317),
            Complex::new(-0.70122531996791548, 1.5685475206682837),
            Complex::new(-1.2959069761255018, -0.56408458111912407),
            Complex::new(0.80142042293998772, -0.15313942732686359),
            Complex::new(0.48871110895462522, -2.5910351676300931),
            Complex::new(2.5341606321914707, -2.0663782886284809),
            Complex::new(0.45027512740123110, -1.6274243463177331),
            Complex::new(-0.57539383024900737, 2.4281348267336553),
            Complex::new(-1.1296731508234275, 0.90398264866704470),
            Complex::new(-0.28563710120446772, -0.47637998868149889),
            Complex::new(-0.58440221723440966, 0.32545120885206402),
            Complex::new(-0.59026260782734252, -0.29890814896817697),
            Complex::new(0.27274978786794912, 2.4996175611869180),
            Complex::new(-0.70330023119372154, -0.95853341627041833),
            Complex::new(0.31510815523064040, -0.27282886655532490),
            Complex::new(-1.0346078437479442, 1.0618637752882218),
            Complex::new(0.65015746875979019, 0.17912169866508931),
            Complex::new(0.89035450285278694, -0.57970274468132166),
            Complex::new(1.6339225000576024E-002, -2.3952402088207974),
            Complex::new(0.33303756158219378, -0.57079592781056654),
            Complex::new(0.22244130892187952, 0.27685893490110408),
            Complex::new(-0.47458924311264428, 1.2988589637548592E-002),
            Complex::new(1.3633048638346255, 0.16949264392950569),
            Complex::new(-0.35756317557586070, 0.57362994646994336),
            Complex::new(-1.0591575896740089, -0.90467331741663981),
            Complex::new(-0.62343016692152187, -0.47254182333647243),
            Complex::new(-4.1375179015731151E-002, -1.4405746313997232),
            Complex::new(-1.1870657309731962, 0.24814596522485438),
            Complex::new(-0.40302903307820248, 1.0136452701079954),
            Complex::new(-0.52516363803653998, 0.56252018882397503),
            Complex::new(0.88670220892379703, -8.8394540644807507E-003),
            Complex::new(0.81042197304080377, -0.30247717946426378),
            Complex::new(1.0014621084444231, 1.3240134188761086),
            Complex::new(-0.40691086281321309, -0.53247309142076349),
            Complex::new(0.99717104313695781, 1.3567383846253696),
            Complex::new(-0.43365445905983463, 1.3229984144889442),
            Complex::new(1.7703807598272381, 0.46455208298342121),
            Complex::new(0.72990142965560356, -1.6605710629343813),
            Complex::new(-0.69601147479901848, -0.98379926699924769),
            Complex::new(-9.0288758923210644E-002, 9.4871429002452989E-002),
            Complex::new(-0.91163509525033404, -0.95184542373586234),
            Complex::new(-1.0277187375192007, 1.7499973452044695),
            Complex::new(-1.7697996686844639, 0.73676561769082727),
            Complex::new(1.3600205778862025, 1.4924824152384053E-002),
        ]
    );

    let mut a = fixtures::complex::matrix_mxn(6, 8);
    let expect = fixtures::complex::matrix_mxn(6, 8);
    complex::gerc(
        6,
        8,
        Complex::new(0.0, 0.0),
        &fixtures::complex::vector(6),
        -1,
        &fixtures::complex::vector(8),
        -1,
        &mut a,
        6,
    );
    capproximately!(a, expect);

    let result = std::panic::catch_unwind(|| {
        complex::gerc(
            6,
            -8,
            Complex::new(0.0, 0.0),
            &vec![],
            1,
            &vec![],
            1,
            &mut vec![],
            6,
        )
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::gerc(
            -6,
            8,
            Complex::new(0.0, 0.0),
            &vec![],
            1,
            &vec![],
            1,
            &mut vec![],
            6,
        )
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::gerc(
            6,
            8,
            Complex::new(0.0, 0.0),
            &vec![],
            1,
            &vec![],
            1,
            &mut vec![],
            2,
        )
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::gerc(
            6,
            8,
            Complex::new(0.0, 0.0),
            &vec![],
            0,
            &vec![],
            1,
            &mut vec![],
            6,
        )
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::gerc(
            6,
            8,
            Complex::new(0.0, 0.0),
            &vec![],
            1,
            &vec![],
            0,
            &mut vec![],
            6,
        )
    });
    assert!(result.is_err());
}

#[test]
fn geru() {
    let mut a = fixtures::complex::matrix_mxn(6, 8);
    let y = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
        Complex::new(-0.9120683669483379, 1.2560188173061),
        Complex::new(-1.4375862408299789, 0.6466743904953449),
    ];

    complex::geru(
        6,
        8,
        Complex::new(0.2, 0.8),
        &fixtures::complex::vector(6),
        1,
        &y,
        1,
        &mut a,
        6,
    );
    capproximately!(
        a,
        vec![
            Complex::new(2.0297612143617068, 0.65171758105023880),
            Complex::new(-0.31722765173412737, -0.22195137259592293),
            Complex::new(0.67838525870385646, 1.3135933588490314),
            Complex::new(0.29660221113116736, -0.36776314720236125),
            Complex::new(0.67085048022165195, 0.96019479554600029),
            Complex::new(-1.5525255484272074, 0.80984860603561071),
            Complex::new(-0.91956134643055376, -0.24522224781145335),
            Complex::new(-0.34337267365494262, -0.84876825882173301),
            Complex::new(3.5884716209125790E-002, -1.3235207357780114),
            Complex::new(2.5134190087092048, -1.2824797785023652),
            Complex::new(0.92423348688558515, -1.4324791531420100),
            Complex::new(-0.85528057562333193, 1.1310756044385850),
            Complex::new(-1.7990710725491952, 0.90733632550094778),
            Complex::new(-0.24780969402183861, -0.38427887944141781),
            Complex::new(0.20431807046274741, 0.36177464242009244),
            Complex::new(0.29590044861983555, -8.3330600887344008E-002),
            Complex::new(-0.12482089968205246, 2.9767875702135211),
            Complex::new(-0.82448777054323985, -0.97825149550145041),
            Complex::new(0.43568331003189087, -5.4877474904060364E-002),
            Complex::new(-1.2375384569168091, 0.25014132261276245),
            Complex::new(-0.22426788508892059, 0.61824327707290649),
            Complex::new(0.37739565968513489, -0.17262350022792816),
            Complex::new(0.13333636522293091, -2.2239003181457520),
            Complex::new(0.80418950319290161, -1.2636144161224365),
            Complex::new(0.19910226600661174, -0.43897941786342354),
            Complex::new(0.66424803995565351, 0.12025746448274954),
            Complex::new(0.70872505335330216, -0.40522610833719930),
            Complex::new(-1.4056841951938006, 0.55246203420080509),
            Complex::new(-1.7086030755539241, -1.5362834602014517),
            Complex::new(0.22396790866494898, 0.41726789091336980),
            Complex::new(-0.24828208812962280, -1.1759958868766818),
            Complex::new(-0.59915959583741274, 0.34047968554477093),
            Complex::new(-0.36587699745669194, 6.5500265815802733E-002),
            Complex::new(-0.49508001953504838, -0.18172001358543421),
            Complex::new(0.90399246920586351, 0.19416080069351449),
            Complex::new(1.0875160462414502, 0.22057782778244783),
            Complex::new(2.1321159056388295, 0.65676437050062852),
            Complex::new(-0.39305997571479284, -8.6017025087555266E-003),
            Complex::new(0.24657229610721754, 1.5163899086758028),
            Complex::new(-1.7869097919659473, 1.2430948361600873),
            Complex::new(2.0642227791536065, -0.83779758149610695),
            Complex::new(0.56860026305061906, -1.1561991079710021),
            Complex::new(0.85772800210885958, -0.88756420944139058),
            Complex::new(-0.96019658037878808, 3.5978134848148047E-003),
            Complex::new(-2.1347229948999118, -0.34176167823825904),
            Complex::new(-2.3793802908792001, 1.6180756075113232),
            Complex::new(-0.66183906951100635, -0.18559064372163170),
            Complex::new(0.96981580269622225, -0.45685234565599903),
        ]
    );

    let mut a = fixtures::complex::matrix_mxn(6, 8);
    complex::geru(
        6,
        8,
        Complex::new(0.2, 0.8),
        &fixtures::complex::vector(6),
        -1,
        &fixtures::complex::vector(8),
        -1,
        &mut a,
        6,
    );
    capproximately!(
        a,
        vec![
            Complex::new(1.0762329818038883, 1.3343172977461006),
            Complex::new(0.57570966920702649, -1.3786971758184945),
            Complex::new(1.6009623633433945E-002, 0.45172643491916009),
            Complex::new(0.30427689638121636, -0.61534079282064869),
            Complex::new(0.28648814967415992, 2.0562213536478549),
            Complex::new(-0.22943803902692839, 0.60174891616896098),
            Complex::new(-0.92071286255973250, -6.9033066835214152E-002),
            Complex::new(1.1599238990566474E-002, -2.0844823115490549),
            Complex::new(-1.5133306773835578, -1.1959050551210406),
            Complex::new(1.4129214686993805, -0.87900000583683158),
            Complex::new(0.80004655638206990, -1.2461504412267883),
            Complex::new(0.34094626223672986, 0.55034719139296251),
            Complex::new(-1.2120527260745289, 0.79528655776482393),
            Complex::new(-0.11221984516143080, -5.2324278013822889E-002),
            Complex::new(-0.14482348710439957, 1.9129169318332417E-002),
            Complex::new(-0.34407749269976445, -0.55961511387532126),
            Complex::new(0.19595214975851932, 2.4159030752180040),
            Complex::new(-0.90449663826714266, -0.54623657520065649),
            Complex::new(0.61292504776955370, 0.12012693498838201),
            Complex::new(-1.6615422282715144, -0.47117343241138221),
            Complex::new(-0.93899822953656631, 1.2865306365243372),
            Complex::new(3.5129813670364829E-004, 0.36279954323693553),
            Complex::new(0.29397641682104048, -2.0925973753786433),
            Complex::new(1.0603985434998502, -2.0613227196364949),
            Complex::new(9.7284846751467785E-002, 0.11172069344247793),
            Complex::new(-0.21112235609010166, 0.65724188116707161),
            Complex::new(2.0311417462261590, -0.29589192058188163),
            Complex::new(1.6457423030846297E-002, 0.17754679974932891),
            Complex::new(-1.1758336062657710, -1.0318578636387787),
            Complex::new(-0.92910096497981542, 0.15384662106908503),
            Complex::new(-0.16827322040571996, -1.6080108071774757),
            Complex::new(-0.91993264527798446, 0.90136418054028800),
            Complex::new(0.27410094386290806, 0.54178477880377995),
            Complex::new(-0.14593846196790139, 0.16092545951632581),
            Complex::new(0.76840262034545570, -0.13779380073408040),
            Complex::new(0.50049769976678371, 0.33262759551857168),
            Complex::new(0.93588906777876768, 1.2374927835462510),
            Complex::new(-0.26887304521127764, -0.19493041442397727),
            Complex::new(1.3470698361170415, 1.1129101360102691),
            Complex::new(-0.23769439837326381, 1.1154791555031165),
            Complex::new(1.7092508665080977, 0.39791647123121499),
            Complex::new(0.56975177917979780, -1.3323882889472414),
            Complex::new(-0.46535950729141939, -0.67946451957474086),
            Complex::new(-0.57583424983129605, -1.0924287446170400),
            Complex::new(-2.1423976801499607, -9.4184032619971658E-002),
            Complex::new(-1.7170046558087777, 2.4799425314411945),
            Complex::new(-1.5547763904521601, 0.97115515950093989),
            Complex::new(1.9233440352540410, -1.1394520623518609),
        ]
    );

    let mut a = fixtures::complex::matrix_mxn(6, 8);
    let expect = fixtures::complex::matrix_mxn(6, 8);
    complex::geru(
        6,
        8,
        Complex::new(0.0, 0.0),
        &fixtures::complex::vector(6),
        -1,
        &fixtures::complex::vector(8),
        -1,
        &mut a,
        6,
    );
    capproximately!(a, expect);

    let result = std::panic::catch_unwind(|| {
        complex::geru(
            6,
            -8,
            Complex::new(0.0, 0.0),
            &vec![],
            1,
            &vec![],
            1,
            &mut vec![],
            6,
        )
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::geru(
            -6,
            8,
            Complex::new(0.0, 0.0),
            &vec![],
            1,
            &vec![],
            1,
            &mut vec![],
            6,
        )
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::geru(
            6,
            8,
            Complex::new(0.0, 0.0),
            &vec![],
            1,
            &vec![],
            1,
            &mut vec![],
            2,
        )
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::geru(
            6,
            8,
            Complex::new(0.0, 0.0),
            &vec![],
            0,
            &vec![],
            1,
            &mut vec![],
            6,
        )
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::geru(
            6,
            8,
            Complex::new(0.0, 0.0),
            &vec![],
            1,
            &vec![],
            0,
            &mut vec![],
            6,
        )
    });
    assert!(result.is_err());
}

#[test]
fn hbmv() {
    let a = matrix::complex::upper_band(
        matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6),
        6,
        6,
        3,
    );
    let x = fixtures::complex::vector(6);
    let mut y = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::hbmv(
        'u',
        6,
        3,
        Complex::new(0.2, 0.8),
        &a,
        4,
        &x,
        1,
        Complex::new(0.3, -0.7),
        &mut y,
        1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(-1.2534610371023731, 9.4703331729354634E-002),
            Complex::new(-0.57417853027717758, -6.9983511992388614E-003),
            Complex::new(-0.24071061122272883, -0.24386226459928881),
            Complex::new(-1.9639075687356289E-003, 0.89891271445795673),
            Complex::new(1.5460009072147165, 0.25069302006213029),
            Complex::new(-1.0003110232032266, -0.49852299275348511)
        ]
    );

    let mut y = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::hbmv(
        'u',
        6,
        3,
        Complex::new(0.2, 0.8),
        &a,
        4,
        &x,
        -1,
        Complex::new(0.0, 0.0),
        &mut y,
        -1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(-1.2736547002392706, -1.2409869471212633),
            Complex::new(-1.5069594023144455, 0.16300769978973884),
            Complex::new(1.3847341260375996, 0.25254847949457943),
            Complex::new(0.41526650594237807, -0.14984859659953415),
            Complex::new(-1.1204385690487695, -0.34592548946596602),
            Complex::new(0.86249408001175165, -0.56909457127832508),
        ]
    );

    let mut y = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];

    complex::hbmv(
        'u',
        6,
        3,
        Complex::new(0.0, 0.0),
        &a,
        4,
        &x,
        -1,
        Complex::new(0.0, 0.0),
        &mut y,
        -1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
        ]
    );

    let mut y = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    let expect = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::hbmv(
        'u',
        6,
        3,
        Complex::new(0.0, 0.0),
        &a,
        4,
        &x,
        -1,
        Complex::new(1.0, 0.0),
        &mut y,
        -1,
    );
    capproximately!(y, expect);

    let mut y = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::hbmv(
        'u',
        6,
        3,
        Complex::new(0.2, 0.8),
        &a,
        4,
        &x,
        -1,
        Complex::new(1.0, 0.0),
        &mut y,
        -1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(-1.9226647624569830, -0.46884477526152213),
            Complex::new(-1.5069594023144455, 0.16300769978973884),
            Complex::new(2.0488698205413836, -0.17226181084202949),
            Complex::new(0.41526650594237807, -0.14984859659953415),
            Complex::new(-0.97666708455456019, 0.65106137653134821),
            Complex::new(0.74474048489799738, -0.84487259666704273),
        ]
    );

    let a = matrix::complex::lower_band(
        matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6),
        6,
        6,
        3,
    );
    let mut y = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::hbmv(
        'l',
        6,
        3,
        Complex::new(0.2, 0.8),
        &a,
        4,
        &x,
        1,
        Complex::new(1.0, 0.0),
        &mut y,
        1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(0.79654686947656517, 0.93230872508886198),
            Complex::new(3.8235731030044517E-002, 3.6497121577365221),
            Complex::new(0.38312451728225538, -3.0675595198095778),
            Complex::new(-1.6321936347511552, -2.1997081054941319),
            Complex::new(3.1857501193367059, 7.7298820669718160E-002),
            Complex::new(0.28503740823625501, 1.0786387593676938),
        ]
    );

    let result = std::panic::catch_unwind(|| {
        complex::hbmv(
            'l',
            -6,
            3,
            Complex::new(0.2, 0.8),
            &a,
            6,
            &x,
            1,
            Complex::new(1.0, 0.0),
            &mut fixtures::complex::vector(6),
            1,
        );
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::hbmv(
            'l',
            6,
            4,
            Complex::new(0.0, 0.0),
            &a,
            4,
            &x,
            1,
            Complex::new(1.0, 0.0),
            &mut fixtures::complex::vector(6),
            1,
        );
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::hbmv(
            'l',
            6,
            -4,
            Complex::new(0.2, 0.8),
            &a,
            6,
            &x,
            1,
            Complex::new(1.0, 0.0),
            &mut fixtures::complex::vector(6),
            1,
        );
    });

    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::hbmv(
            'l',
            6,
            4,
            Complex::new(0.2, 0.8),
            &a,
            5,
            &x,
            0,
            Complex::new(1.0, 0.0),
            &mut fixtures::complex::vector(6),
            1,
        );
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::hbmv(
            'l',
            6,
            4,
            Complex::new(0.0, 0.0),
            &a,
            5,
            &x,
            1,
            Complex::new(1.0, 0.0),
            &mut fixtures::complex::vector(6),
            0,
        );
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::hbmv(
            'x',
            6,
            4,
            Complex::new(0.0, 0.0),
            &a,
            5,
            &x,
            1,
            Complex::new(1.0, 0.0),
            &mut fixtures::complex::vector(6),
            1,
        );
    });
    assert!(result.is_err());
}

#[test]
fn hemv() {
    let mut a = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut a, 6, 6, 0.0);
    let x = fixtures::complex::vector(6);
    let mut y = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::hemv(
        'u',
        6,
        Complex::new(0.2, 0.8),
        &a,
        6,
        &x,
        1,
        Complex::new(0.3, -0.7),
        &mut y,
        1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(-1.1978771580109713, 0.19256169948001889),
            Complex::new(-0.64406391232734050, 0.12304172454730536),
            Complex::new(-0.24071061122272883, -0.24386226459928881),
            Complex::new(-1.9639075687356289E-003, 0.89891271445795673),
            Complex::new(1.5927182864486373, 0.26326773843256362),
            Complex::new(-1.0745919479739299, -0.10110828501746819),
        ]
    );

    let mut y = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::hemv(
        'u',
        6,
        Complex::new(-0.12, 0.88),
        &a,
        6,
        &x,
        -1,
        Complex::new(-0.43, 0.57),
        &mut y,
        -1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(-0.53727175687523987, -2.3112656523209361),
            Complex::new(-1.5860492969147173, -0.43291167355929333),
            Complex::new(1.2402782186527825, 1.3676678793636570),
            Complex::new(0.47520593481607637, 1.6258004364783057E-002),
            Complex::new(-1.1206760713901169, -0.92993994108923006),
            Complex::new(1.6699524223689810, -0.23413210623927322)
        ]
    );

    let mut a = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_upper(&mut a, 6, 6, 0.0);
    let mut y = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::hemv(
        'l',
        6,
        Complex::new(-0.12, 0.88),
        &a,
        6,
        &x,
        -1,
        Complex::new(-0.43, 0.57),
        &mut y,
        -1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(-0.82109209135990002, -2.8288933164643106),
            Complex::new(1.8186627440037402, -1.9201765078819051),
            Complex::new(-4.8025707150157597, 2.6701570029481525),
            Complex::new(2.3815496358250265, 0.20742891107446010),
            Complex::new(-1.7405950102646242, 2.9574110901705071),
            Complex::new(3.6227184273789095, 3.3816196677830903),
        ]
    );

    let mut y = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];

    complex::hemv(
        'l',
        6,
        Complex::new(0.0, 0.0),
        &a,
        6,
        &x,
        -1,
        Complex::new(0.0, 0.0),
        &mut y,
        -1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
        ]
    );

    let mut y = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];

    complex::hemv(
        'l',
        6,
        Complex::new(0.0, 0.0),
        &a,
        6,
        &x,
        -1,
        Complex::new(1.0, 0.0),
        &mut y,
        -1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(-0.6490100777088978, 0.7721421858045301),
            Complex::new(0.0, 0.0),
            Complex::new(0.6641356998941105, -0.4248102833772871),
            Complex::new(0.0, 0.0),
            Complex::new(0.14377148075806995, 0.9969868609091059),
            Complex::new(-0.11775359816595128, -0.27577802908802723),
        ]
    );

    let mut y = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::hemv(
        'l',
        6,
        Complex::new(0.0, 1.0),
        &a,
        6,
        &x,
        -1,
        Complex::new(1.0, 0.0),
        &mut y,
        -1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(-1.7089394711323604, -1.5002952495366637),
            Complex::new(1.7368180095711416, -2.4188575862567556),
            Complex::new(-4.3244119692369942, 2.6519589765813314),
            Complex::new(2.6884573545075430, -0.13089314190862256),
            Complex::new(-0.59245062296078643, 4.8521145635536458),
            Complex::new(4.1985782666268578, 2.9198979594665389),
        ]
    );

    let result = std::panic::catch_unwind(|| {
        complex::hemv(
            'l',
            -6,
            Complex::new(0.2, 0.8),
            &a,
            6,
            &x,
            1,
            Complex::new(1.0, 0.0),
            &mut fixtures::complex::vector(6),
            1,
        );
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::hemv(
            'l',
            6,
            Complex::new(0.0, 0.0),
            &a,
            5,
            &x,
            1,
            Complex::new(1.0, 0.0),
            &mut fixtures::complex::vector(6),
            1,
        );
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::hemv(
            'l',
            6,
            Complex::new(0.0, 0.0),
            &a,
            6,
            &x,
            0,
            Complex::new(1.0, 0.0),
            &mut fixtures::complex::vector(6),
            1,
        );
    });

    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::hemv(
            'l',
            6,
            Complex::new(0.2, 0.8),
            &a,
            5,
            &x,
            1,
            Complex::new(1.0, 0.0),
            &mut fixtures::complex::vector(6),
            0,
        );
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::hemv(
            'x',
            6,
            Complex::new(0.0, 0.0),
            &a,
            5,
            &x,
            1,
            Complex::new(1.0, 0.0),
            &mut fixtures::complex::vector(6),
            1,
        );
    });
    assert!(result.is_err());
}

#[test]
fn her() {
    let mut a = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut a, 6, 6, 0.0);
;
    let x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];

    complex::her('u', 6, 0.2, &x, 1, &mut a, 6);
    capproximately!(
        a,
        vec![
            Complex::new(1.4664377569938685, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.92856705188751221, -0.83204329013824463),
            Complex::new(-0.29472044110298157, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-1.2994659767673966, -0.17990848227932066),
            Complex::new(-0.28946158289909363, 0.26613736152648926),
            Complex::new(-0.17490710560589112, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.43568331003189087, 0.25014132261276245),
            Complex::new(-1.2375384569168091, 0.61824327707290649),
            Complex::new(-0.22426788508892059, -0.17262350022792816),
            Complex::new(0.37739565968513489, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(7.8194520501826237E-002, 0.14056783081885876),
            Complex::new(0.50360798835754395, -0.94064915180206299),
            Complex::new(1.0201601128637381, -0.26046736155203748),
            Complex::new(-0.69095385074615479, -0.81496870517730713),
            Complex::new(-1.0816686910303539, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.26300986834671181, 0.31196009025279664),
            Complex::new(-0.54288828372955322, 0.24841265380382538),
            Complex::new(-0.42552053432548659, 0.11192357318741848),
            Complex::new(-0.64947164058685303, 1.9156390801072121E-002),
            Complex::new(0.66837539491488840, 0.24178842029114961),
            Complex::new(1.1698956214925516, 0.0000000000000000),
        ]
    );

    let mut a = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut a, 6, 6, 0.0);
    complex::her('u', 6, 0.2, &x, -1, &mut a, 6);
    capproximately!(
        a,
        vec![
            Complex::new(1.2809381210347879, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.98694238844082449, -0.81649333557626003),
            Complex::new(-9.1789827934116691E-002, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-1.1476570367813110, -0.22732868790626526),
            Complex::new(-0.28946158289909363, 0.26613736152648926),
            Complex::new(-0.29921510815620422, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.44347310562049730, 0.20350592804806661),
            Complex::new(-1.3031477589548044, 0.76288531337431831),
            Complex::new(-0.22426788508892059, -0.17262350022792816),
            Complex::new(0.50170366223544804, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-5.7106774300336838E-002, -1.1045478284358978E-002),
            Complex::new(0.50360798835754395, -0.94064915180206299),
            Complex::new(1.0857694149017334, -0.11582532525062561),
            Complex::new(-0.69095385074615479, -0.81496870517730713),
            Complex::new(-1.2845993041992188, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.26300986834671181, 0.41992218389805175),
            Complex::new(-0.40758698892739015, 9.6799344700607637E-002),
            Complex::new(-0.43331032991409302, 6.5288178622722626E-002),
            Complex::new(-0.80128058057293872, -2.8263814825872494E-002),
            Complex::new(0.72675073146820068, 0.25733837485313416),
            Complex::new(1.3553952574516321, 0.0000000000000000),
        ]
    );

    let mut a = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut a, 6, 6, 0.0);
    complex::her('u', 6, 0.0, &x, -1, &mut a, 6);
    capproximately!(
        a,
        vec![
            Complex::new(1.2629542350769043, -0.42951309680938721),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.92856705188751221, -0.83204329013824463),
            Complex::new(-0.29472044110298157, -1.1665705442428589),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-1.1476570367813110, -0.22732868790626526),
            Complex::new(-0.28946158289909363, 0.26613736152648926),
            Complex::new(-0.29921510815620422, -0.37670272588729858),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.43568331003189087, 0.25014132261276245),
            Complex::new(-1.2375384569168091, 0.61824327707290649),
            Complex::new(-0.22426788508892059, -0.17262350022792816),
            Complex::new(0.37739565968513489, -2.2239003181457520),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-5.7106774300336838E-002, -1.1045478284358978E-002),
            Complex::new(0.50360798835754395, -0.94064915180206299),
            Complex::new(1.0857694149017334, -0.11582532525062561),
            Complex::new(-0.69095385074615479, -0.81496870517730713),
            Complex::new(-1.2845993041992188, 0.24226348102092743),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.23570655286312103, 0.36594113707542419),
            Complex::new(-0.54288828372955322, 0.24841265380382538),
            Complex::new(-0.43331032991409302, 6.5288178622722626E-002),
            Complex::new(-0.64947164058685303, 1.9156390801072121E-002),
            Complex::new(0.72675073146820068, 0.25733837485313416),
            Complex::new(1.1519117355346680, 0.0000000000000000),
        ]
    );

    let mut a = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_upper(&mut a, 6, 6, 0.0);
    complex::her('l', 6, 1.234, &x, -1, &mut a, 6);
    capproximately!(
        a,
        vec![
            Complex::new(1.3739148068679738, 0.0000000000000000),
            Complex::new(-0.68640916889390469, 1.1423609224868452),
            Complex::new(1.3297992944717407, -0.27934628725051880),
            Complex::new(1.3204923838408580, 2.0456434716758203),
            Complex::new(0.41464143991470337, 0.56074607372283936),
            Complex::new(-1.7084114627576494, -0.78584701720584593),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.95736139059139025, 0.0000000000000000),
            Complex::new(-5.7671726681292057E-003, -1.0655906200408936),
            Complex::new(1.9998439338703389, -2.4562234231403806),
            Complex::new(0.76359343528747559, 1.1565370559692383),
            Complex::new(3.5799691038578052E-002, 1.7675011834827594),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.29921510815620422, 0.0000000000000000),
            Complex::new(-0.41151082515716553, 2.4413645267486572),
            Complex::new(0.25222346186637878, -0.79533910751342773),
            Complex::new(-0.89192110300064087, -5.4877474904060364E-002),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(1.1443760038382795, 0.0000000000000000),
            Complex::new(0.13333636522293091, -1.2636144161224365),
            Complex::new(-0.13247161795193985, 0.65131154232107824),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-1.2845993041992188, 0.0000000000000000),
            Complex::new(4.6726170927286148E-002, -1.4250984191894531),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(2.4074050140643375, 0.0000000000000000),
        ]
    );

    let result = std::panic::catch_unwind(|| {
        complex::her('l', -6, 5.0, &x, 1, &mut vec![], 6);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::her('l', 6, 1.2, &x, 1, &mut vec![], 5);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::her('l', 6, 0.2, &x, 0, &mut vec![], 6);
    });

    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::her('x', 6, 0.1, &x, 1, &mut vec![], 6);
    });

    assert!(result.is_err());
}

#[test]
fn her2() {
    let mut a = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut a, 6, 6, 0.0);
    let y = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    let x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];

    complex::her2('u', 6, Complex::new(0.2, 0.8), &x, 1, &y, 1, &mut a, 6);
    capproximately!(
        a,
        vec![
            Complex::new(1.6699212789108326, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.75941169595474800, -0.95266433991192034),
            Complex::new(-0.29472044110298157, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-1.4512749167534824, -0.13248827665237595),
            Complex::new(-0.44376975470840946, 0.21556829407581407),
            Complex::new(-5.0599103055578076E-002, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.43568331003189087, 0.25014132261276245),
            Complex::new(-1.2375384569168091, 0.61824327707290649),
            Complex::new(-0.22426788508892059, -0.17262350022792816),
            Complex::new(0.37739565968513489, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.21349581530398931, 0.29218113992207650),
            Complex::new(0.52621022235371806, -0.73440941757238531),
            Complex::new(0.95455081082574245, -0.40510939785344929),
            Complex::new(-0.69095385074615479, -0.81496870517730713),
            Complex::new(-0.87873807786148905, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.29031318383030258, 0.25797904343016909),
            Complex::new(-0.53358655517175724, 0.18735347000462102),
            Complex::new(-0.41773073873688016, 0.15855896775211434),
            Complex::new(-0.64947164058685303, 1.9156390801072121E-002),
            Complex::new(0.61000005836157611, 0.22623846572916503),
            Complex::new(1.1878795074504351, 0.0000000000000000),
        ]
    );
    let mut a = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_upper(&mut a, 6, 6, 0.0);
    complex::her2('l', 6, Complex::new(0.2, 0.8), &x, -1, &y, -1, &mut a, 6);
    capproximately!(
        a,
        vec![
            Complex::new(1.2989220069926715, 0.0000000000000000),
            Complex::new(-0.44298403029771038, 1.2072042290596245),
            Complex::new(1.3297992944717407, -0.27934628725051880),
            Complex::new(1.2880089382154818, 1.8511738881894504),
            Complex::new(0.42394316847249935, 0.49968688992363497),
            Complex::new(-1.5945566441278871, -0.56074606567017271),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.11114078523474813, 0.0000000000000000),
            Complex::new(-5.7671726681292057E-003, -1.0655906200408936),
            Complex::new(2.2734347066997658, -1.8530661685119423),
            Complex::new(0.78619566928364970, 1.3627767901989158),
            Complex::new(-0.52840667391114637, 1.1352737230419456),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.29921510815620422, 0.0000000000000000),
            Complex::new(-0.41151082515716553, 2.4413645267486572),
            Complex::new(0.25222346186637878, -0.79533910751342773),
            Complex::new(-0.89192110300064087, -5.4877474904060364E-002),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.62601166478576098, 0.0000000000000000),
            Complex::new(-2.0971806586384922E-002, -1.3141834835731117),
            Complex::new(0.50057162322073023, 0.45356929690452408),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-1.2845993041992188, 0.0000000000000000),
            Complex::new(0.21588152686005035, -1.5457194689631288),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(1.5588787793685963, 0.0000000000000000),
        ]
    );

    let mut a = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut a, 6, 6, 0.0);
    let mut expect = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut expect, 6, 6, 0.0);
    complex::her2('l', 6, Complex::new(0.0, 0.0), &x, -1, &y, -1, &mut a, 6);
    capproximately!(a, expect);

    let result = std::panic::catch_unwind(|| {
        complex::her2(
            'x',
            6,
            Complex::new(0.0, 0.0),
            &x,
            -1,
            &y,
            -1,
            &mut vec![],
            6,
        );
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::her2(
            'l',
            -6,
            Complex::new(0.0, 0.0),
            &x,
            -1,
            &y,
            -1,
            &mut vec![],
            6,
        );
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::her2(
            'l',
            6,
            Complex::new(0.0, 0.0),
            &x,
            0,
            &y,
            -1,
            &mut vec![],
            6,
        );
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::her2(
            'l',
            6,
            Complex::new(0.0, 0.0),
            &x,
            -1,
            &y,
            0,
            &mut vec![],
            6,
        );
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::her2(
            'x',
            6,
            Complex::new(0.0, 0.0),
            &x,
            -1,
            &y,
            -1,
            &mut vec![],
            4,
        );
    });
    assert!(result.is_err());
}

#[test]
fn hpmv() {
    let mut ap = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut ap, 6, 6, 0.0);
    let ap = matrix::complex::pack_upper(ap, 6, 6, 5);
    let mut y = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];

    let x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];

    complex::hpmv(
        'u',
        6,
        Complex::new(0.2, -0.8),
        &ap,
        &x,
        1,
        Complex::new(0.3, -0.7),
        &mut y,
        1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(1.0565062967743115, 2.1785481727394473),
            Complex::new(-0.36064336254715390, -0.89535796823641833),
            Complex::new(0.26097073349390854, -1.1070832376210971),
            Complex::new(0.18952760417765735, -0.48786643064345248),
            Complex::new(-0.69670093526894639, -0.67058224989704085),
            Complex::new(0.29360012462781360, -0.16342854670718146),
        ]
    );

    let mut y = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];

    complex::hpmv(
        'l',
        6,
        Complex::new(0.2, -0.8),
        &ap,
        &x,
        -1,
        Complex::new(0.3, -0.7),
        &mut y,
        -1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(1.0792789987189280, 1.7817909884398611),
            Complex::new(1.1638084092421319, -1.2174375200025294),
            Complex::new(-9.8419729711608492E-002, -0.92170300747050660),
            Complex::new(-1.5742564828419501, -0.99686730476185714),
            Complex::new(1.4777738731392653, 0.47151670288681441),
            Complex::new(-1.3963398509841487, 1.3006511215030607),
        ]
    );

    let mut ap = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut ap, 6, 6, 0.0);
    let ap = matrix::complex::pack_upper(ap, 6, 6, 5);
    let mut y = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];

    complex::hpmv(
        'l',
        6,
        Complex::new(0.0, 0.0),
        &ap,
        &x,
        -1,
        Complex::new(1.0, 0.0),
        &mut y,
        -1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(-0.6490100777088978, 0.7721421858045301),
            Complex::new(-0.11916876241803812, -0.21951562675343952),
            Complex::new(0.6641356998941105, -0.4248102833772871),
            Complex::new(0.0, 0.0),
            Complex::new(0.14377148075806995, 0.9969868609091059),
            Complex::new(-0.11775359816595128, -0.27577802908802723),
        ]
    );

    let mut y = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];

    complex::hpmv(
        'l',
        6,
        Complex::new(0.0, 0.0),
        &ap,
        &x,
        -1,
        Complex::new(0.0, 0.0),
        &mut y,
        -1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0)
        ]
    );

    let mut ap = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut ap, 6, 6, 0.0);
    let ap = matrix::complex::pack_upper(ap, 6, 6, 5);
    let mut y = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];

    complex::hpmv(
        'l',
        6,
        Complex::new(1.0, 0.0),
        &ap,
        &x,
        -1,
        Complex::new(1.0, 0.0),
        &mut y,
        -1,
    );
    capproximately!(
        y,
        vec![
            Complex::new(-1.7225049442869467, 1.9573690212447206),
            Complex::new(1.7317793933958838, 1.0092722795253204),
            Complex::new(1.0515376006342987, -0.52202732662423190),
            Complex::new(0.70976844124876681, -2.1452626845418177),
            Complex::new(3.9215297000845251E-002, 1.9440654388525900),
            Complex::new(-1.9918115746842791, -1.2672249544431136)
        ]
    );

    let result = std::panic::catch_unwind(|| {
        complex::hpmv(
            'x',
            6,
            Complex::new(0.0, 0.0),
            &vec![],
            &vec![],
            -1,
            Complex::new(0.0, 0.0),
            &mut vec![],
            -1,
        );
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::hpmv(
            'l',
            -6,
            Complex::new(0.0, 0.0),
            &vec![],
            &vec![],
            -1,
            Complex::new(0.0, 0.0),
            &mut vec![],
            -1,
        );
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::hpmv(
            'l',
            6,
            Complex::new(0.0, 0.0),
            &vec![],
            &vec![],
            0,
            Complex::new(0.0, 0.0),
            &mut vec![],
            -1,
        );
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::hpmv(
            'l',
            6,
            Complex::new(0.0, 0.0),
            &vec![],
            &vec![],
            -1,
            Complex::new(0.0, 0.0),
            &mut vec![],
            0,
        );
    });
    assert!(result.is_err());
}

#[test]
fn hpr() {
    let mut ap = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut ap, 6, 6, 0.0);
    let mut ap = matrix::complex::pack_upper(ap, 6, 6, 5);
    let x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::hpr('u', 6, 0.2, &x, 1, &mut ap);
    capproximately!(
        ap,
        vec![
            Complex::new(1.4664377569938685, 0.0000000000000000),
            Complex::new(-0.92856705188751221, -0.83204329013824463),
            Complex::new(-0.29472044110298157, 0.0000000000000000),
            Complex::new(-1.2994659767673966, -0.17990848227932066),
            Complex::new(-0.28946158289909363, 0.26613736152648926),
            Complex::new(-0.17490710560589112, 0.0000000000000000),
            Complex::new(0.43568331003189087, 0.25014132261276245),
            Complex::new(-1.2375384569168091, 0.61824327707290649),
            Complex::new(-0.22426788508892059, -0.17262350022792816),
            Complex::new(0.37739565968513489, 0.0000000000000000),
            Complex::new(7.8194520501826237E-002, 0.14056783081885876),
            Complex::new(0.50360798835754395, -0.94064915180206299),
            Complex::new(1.0201601128637381, -0.26046736155203748),
            Complex::new(-0.69095385074615479, -0.81496870517730713),
            Complex::new(-1.0816686910303539, 0.0000000000000000),
            Complex::new(-0.26300986834671181, 0.31196009025279664),
            Complex::new(-0.54288828372955322, 0.24841265380382538),
            Complex::new(-0.42552053432548659, 0.11192357318741848),
            Complex::new(-0.64947164058685303, 1.9156390801072121E-002),
            Complex::new(0.66837539491488840, 0.24178842029114961),
            Complex::new(1.1698956214925516, 0.0000000000000000),
        ]
    );
    let mut ap = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut ap, 6, 6, 0.0);
    let mut ap = matrix::complex::pack_upper(ap, 6, 6, 5);
    complex::hpr('l', 6, 0.2, &x, -1, &mut ap);
    capproximately!(
        ap,
        vec![
            Complex::new(1.2809381210347879, 0.0000000000000000),
            Complex::new(-0.98694238844082449, -0.84759324470022923),
            Complex::new(-0.29472044110298157, -1.1665705442428589),
            Complex::new(-1.1398672411927047, -0.18069329334156942),
            Complex::new(-0.28946158289909363, 0.26613736152648926),
            Complex::new(-0.32651842363979500, -0.43068377270992614),
            Complex::new(0.63861392320075572, 0.0000000000000000),
            Complex::new(-1.2375384569168091, 0.61824327707290649),
            Complex::new(-0.28987718712691601, -0.31726553652934003),
            Complex::new(0.37739565968513489, -2.2239003181457520),
            Complex::new(7.8194520501826237E-002, 0.14056783081885876),
            Complex::new(0.50360798835754395, 0.0000000000000000),
            Complex::new(1.0857694149017334, -0.11582532525062561),
            Complex::new(-0.69095385074615479, -0.81496870517730713),
            Complex::new(-1.2845993041992188, 0.24226348102092743),
            Complex::new(-0.11139855031280793, 0.0000000000000000),
            Complex::new(-0.54288828372955322, 0.24841265380382538),
            Complex::new(-0.58511926990017871, 0.11270838424966724),
            Complex::new(-0.64947164058685303, 0.0000000000000000),
            Complex::new(0.72675073146820068, 0.25733837485313416),
            Complex::new(1.3553952574516321, 0.0000000000000000),
        ]
    );
    let mut ap = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut ap, 6, 6, 0.0);
    let mut ap = matrix::complex::pack_upper(ap, 6, 6, 5);

    let mut expect = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut expect, 6, 6, 0.0);
    let expect = matrix::complex::pack_upper(expect, 6, 6, 5);

    complex::hpr('l', 6, 0.0, &x, -1, &mut ap);
    capproximately!(ap, expect);

    let result = std::panic::catch_unwind(|| {
        complex::hpr('x', 6, 0.0, &x, 1, &mut vec![]);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::hpr('l', -6, 0.0, &x, 1, &mut vec![]);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::hpr('l', 6, 0.0, &x, 0, &mut vec![]);
    });
    assert!(result.is_err());
}

#[test]
fn hpr2() {
    let mut ap = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut ap, 6, 6, 0.0);
    let mut ap = matrix::complex::pack_upper(ap, 6, 6, 5);
    let y = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    let x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::hpr2('u', 6, Complex::new(0.2, 0.8), &x, 1, &y, 1, &mut ap);
    capproximately!(
        ap,
        vec![
            Complex::new(1.6699212789108326, 0.0000000000000000),
            Complex::new(-0.75941169595474800, -0.95266433991192034),
            Complex::new(-0.29472044110298157, 0.0000000000000000),
            Complex::new(-1.4512749167534824, -0.13248827665237595),
            Complex::new(-0.44376975470840946, 0.21556829407581407),
            Complex::new(-5.0599103055578076E-002, 0.0000000000000000),
            Complex::new(0.43568331003189087, 0.25014132261276245),
            Complex::new(-1.2375384569168091, 0.61824327707290649),
            Complex::new(-0.22426788508892059, -0.17262350022792816),
            Complex::new(0.37739565968513489, 0.0000000000000000),
            Complex::new(0.21349581530398931, 0.29218113992207650),
            Complex::new(0.52621022235371806, -0.73440941757238531),
            Complex::new(0.95455081082574245, -0.40510939785344929),
            Complex::new(-0.69095385074615479, -0.81496870517730713),
            Complex::new(-0.87873807786148905, 0.0000000000000000),
            Complex::new(-0.29031318383030258, 0.25797904343016909),
            Complex::new(-0.53358655517175724, 0.18735347000462102),
            Complex::new(-0.41773073873688016, 0.15855896775211434),
            Complex::new(-0.64947164058685303, 1.9156390801072121E-002),
            Complex::new(0.61000005836157611, 0.22623846572916503),
            Complex::new(1.1878795074504351, 0.0000000000000000),
        ]
    );
    let mut ap = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_upper(&mut ap, 6, 6, 0.0);
    let mut ap = matrix::complex::pack_lower(ap, 6, 6, 5);
    complex::hpr2('l', 6, Complex::new(0.2, 0.8), &x, -1, &y, -1, &mut ap);
    capproximately!(
        ap,
        vec![
            Complex::new(1.2989220069926715, 0.0000000000000000),
            Complex::new(-0.44298403029771038, 1.2072042290596245),
            Complex::new(1.3297992944717407, -0.27934628725051880),
            Complex::new(1.2880089382154818, 1.8511738881894504),
            Complex::new(0.42394316847249935, 0.49968688992363497),
            Complex::new(-1.5945566441278871, -0.56074606567017271),
            Complex::new(0.11114078523474813, 0.0000000000000000),
            Complex::new(-5.7671726681292057E-003, -1.0655906200408936),
            Complex::new(2.2734347066997658, -1.8530661685119423),
            Complex::new(0.78619566928364970, 1.3627767901989158),
            Complex::new(-0.52840667391114637, 1.1352737230419456),
            Complex::new(-0.29921510815620422, 0.0000000000000000),
            Complex::new(-0.41151082515716553, 2.4413645267486572),
            Complex::new(0.25222346186637878, -0.79533910751342773),
            Complex::new(-0.89192110300064087, -5.4877474904060364E-002),
            Complex::new(0.62601166478576098, 0.0000000000000000),
            Complex::new(-2.0971806586384922E-002, -1.3141834835731117),
            Complex::new(0.50057162322073023, 0.45356929690452408),
            Complex::new(-1.2845993041992188, 0.0000000000000000),
            Complex::new(0.21588152686005035, -1.5457194689631288),
            Complex::new(1.5588787793685963, 0.0000000000000000),
        ]
    );
    let mut ap = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_upper(&mut ap, 6, 6, 0.0);
    let mut ap = matrix::complex::pack_lower(ap, 6, 6, 5);
    let mut expect = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_upper(&mut expect, 6, 6, 0.0);
    let expect = matrix::complex::pack_lower(expect, 6, 6, 5);

    complex::hpr2('l', 6, Complex::new(0.0, 0.0), &x, -1, &y, -1, &mut ap);
    capproximately!(ap, expect);

    let result = std::panic::catch_unwind(|| {
        complex::hpr2('x', 6, Complex::new(0.0, 0.0), &x, 1, &y, 1, &mut vec![]);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::hpr2('l', -6, Complex::new(0.0, 0.0), &x, 1, &y, 1, &mut vec![]);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::hpr2('l', 6, Complex::new(0.0, 0.0), &x, 0, &y, 1, &mut vec![]);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::hpr2('l', 6, Complex::new(0.0, 0.0), &x, 1, &y, 0, &mut vec![]);
    });
    assert!(result.is_err());
}

#[test]
fn tbmv() {
    let a = fixtures::complex::matrix_mxn(6, 6);
    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];

    complex::tbmv('u', 'n', 'n', 6, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.89893773605253458, 1.1357840738783445),
            Complex::new(-0.45151985700522917, 0.29850125888377299),
            Complex::new(-0.71516509118603144, 0.65794798179721126),
            Complex::new(1.2587631018024554, 1.0538850956187753),
            Complex::new(0.13566746748087732, -0.61527830905039860),
            Complex::new(9.4482665585330361E-002, 0.17142208883575871),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tbmv('u', 'n', 'u', 6, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.93782339751986754, 0.74413133490913808),
            Complex::new(-0.45151985700522917, 0.29850125888377299),
            Complex::new(0.38229682533410747, 0.30850538482755907),
            Complex::new(1.2587631018024554, 1.0538850956187753),
            Complex::new(0.26330208478911077, 1.0872328501044080),
            Complex::new(-0.11775359511375427, -0.27577802538871765),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tbmv('l', 'n', 'n', 6, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-1.5857588772442313, 0.33125815615156995),
            Complex::new(0.54337390686068332, 2.6859788728962997E-002),
            Complex::new(-2.2279378005181627, 1.2632509704356920),
            Complex::new(-0.89893773605253458, 1.1357840738783445),
            Complex::new(-0.45151985700522917, 0.29850125888377299),
            Complex::new(-0.71516509118603144, 0.65794798179721126),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tbmv('l', 'n', 'u', 6, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.64901006221771240, 0.77214217185974121),
            Complex::new(0.54337390686068332, 2.6859788728962997E-002),
            Complex::new(-1.1550642750183755, -0.20168802073263237),
            Complex::new(-0.89893773605253458, 1.1357840738783445),
            Complex::new(5.8109940652350645E-002, 1.3008478443970510),
            Complex::new(-0.46766315226263799, 0.14935680643823457),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tbmv('u', 't', 'n', 6, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.61012440075037944, 1.1637949108289476),
            Complex::new(0.90450126675144427, 0.75266294427094904),
            Complex::new(-6.9932675283407519E-002, -0.15133458277470524),
            Complex::new(-0.12669784107884396, 0.87789422422181218),
            Complex::new(0.34590867545458170, -0.92679784529900511),
            Complex::new(-0.97741676538200184, -1.0711961221936419),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tbmv('u', 't', 'u', 6, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.64901006221771240, 0.77214217185974121),
            Complex::new(0.90450126675144427, 0.75266294427094904),
            Complex::new(1.0275292412367314, -0.50077717974435743),
            Complex::new(-0.12669784107884396, 0.87789422422181218),
            Complex::new(0.47354329276281515, 0.77571331385580145),
            Complex::new(-1.1896530260810865, -1.5183962364181183),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tbmv('l', 't', 'n', 6, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.17654735879280281, 0.58874771059708308),
            Complex::new(0.85890551240220603, 1.8168168506139182),
            Complex::new(-0.77252134008506745, 0.93792187448835485),
            Complex::new(-0.23040409445892918, -1.2087986124495775),
            Complex::new(-0.42820606451121151, -0.14294309133745720),
            Complex::new(-0.36525553403714772, 0.23281314997025904),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tbmv('l', 't', 'u', 6, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.76020145623371604, 1.0296317263052543),
            Complex::new(0.85890551240220603, 1.8168168506139182),
            Complex::new(0.30035218541471975, -0.52701711667996953),
            Complex::new(-0.23040409445892918, -1.2087986124495775),
            Complex::new(8.1423733146368305E-002, 0.85940349417582085),
            Complex::new(-0.11775359511375427, -0.27577802538871765),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tbmv('u', 'c', 'n', 6, 3, &a, 6, &mut x, -1);
    capproximately!(
        x,
        vec![
            Complex::new(-4.4081662135385002E-002, -0.47152180004162414),
            Complex::new(1.4701333254503495, 5.4962968975589277E-002),
            Complex::new(0.35926976559701229, -1.4420628114452261),
            Complex::new(0.31909425240906009, -0.28351812300138768),
            Complex::new(-0.39426565004330849, 2.4148318973788592),
            Complex::new(-7.2795562644385470E-002, -0.38380208237829727),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tbmv('u', 'c', 'u', 6, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.64901006221771240, 0.77214217185974121),
            Complex::new(-0.89701536056711850, -0.76156909872986789),
            Complex::new(0.67646910762479351, -0.79585439150818971),
            Complex::new(-0.73671647846003907, -1.4532698964078605E-002),
            Complex::new(0.48292775843658564, 0.79038470663881988),
            Complex::new(0.51646822672704484, 0.30309190830210242),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tbmv('l', 'c', 'n', 6, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.30354170139129621, 0.23178789271847933),
            Complex::new(-0.55893308970669509, 3.2283992379330400),
            Complex::new(-0.74100210864921268, -0.33250461678806609),
            Complex::new(-7.2624688204718790E-002, -1.1351242539598538),
            Complex::new(0.29318211076532813, -0.24869434973481597),
            Complex::new(0.42076612202015307, -0.10280777453071321),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tbmv('l', 'c', 'u', 6, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.29188722469987738, -0.61517223133368759),
            Complex::new(-0.55893308970669509, 3.2283992379330400),
            Complex::new(1.0387957612083905, -0.69225924429680141),
            Complex::new(-7.2624688204718790E-002, -1.1351242539598538),
            Complex::new(8.7515933527837264E-002, 0.85680220462035250),
            Complex::new(-0.11775359511375427, -0.27577802538871765),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tbmv('l', 'c', 'u', 0, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.64901006221771240, 0.77214217185974121),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.66413569450378418, -0.42481029033660889),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.14377148449420929, 0.99698686599731445),
            Complex::new(-0.11775359511375427, -0.27577802538871765),
        ]
    );

    let result = std::panic::catch_unwind(|| {
        complex::tbmv('x', 't', 'u', 6, 5, &a, 6, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::tbmv('u', 'x', 'u', 6, 5, &a, 6, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::tbmv('u', 't', 'c', 6, 5, &a, 6, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::tbmv('u', 't', 'n', -6, 5, &a, 6, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::tbmv('u', 't', 'n', 6, -5, &a, 6, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::tbmv('u', 't', 'n', 6, 5, &a, 4, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::tbmv('u', 't', 'n', 6, 5, &a, 4, &mut vec![], 0);
    });
    assert!(result.is_err());
}

#[test]
fn tbsv() {
    let a = fixtures::complex::matrix_mxn(6, 6);
    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(0.0, 0.0),
    ];

    complex::tbsv('u', 'n', 'n', 6, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-8.7140181184168313, 0.53974179877620077),
            Complex::new(1.8661327422332243, 3.9027499877685838),
            Complex::new(-3.6580339360462011, 8.8097913188003023),
            Complex::new(5.6188869330661131, 4.3092860094409744),
            Complex::new(-0.43765579556553402, -1.3695491241654123),
            Complex::new(0.0000000000000000, 0.0000000000000000),
        ]
    );
    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tbsv('u', 'n', 'u', 6, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(3.3386606552999880, -0.30480944658021930),
            Complex::new(-0.45872546160875621, -2.5016404952691627),
            Complex::new(0.10040131610044789, -1.4827193985502860),
            Complex::new(5.6878452703088889E-002, -1.4873152158496024),
            Complex::new(2.4240884199307811E-002, 0.90674088189022095),
            Complex::new(-0.11775359511375427, -0.27577802538871765),
        ]
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tbsv('l', 'n', 'n', 6, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.19681361062735658, 0.14043312545782299),
            Complex::new(-0.43666415573931500, -0.12517747762468392),
            Complex::new(1.9908862356588675, -0.49306176041025807),
            Complex::new(-1.1707260559995081, -4.5452002398639015),
            Complex::new(-0.36248447218257690, 0.31276418173763804),
        ]
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.0, 0.0),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];

    complex::tbsv('l', 'n', 'u', 6, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.11916876584291458, -0.21951562166213989),
            Complex::new(0.14752502884958729, -0.16384931285891380),
            Complex::new(1.4363127452141604, -0.57315625854518182),
            Complex::new(2.2989023266518371, 0.24098851995891635),
            Complex::new(-1.1879615574537852, -1.4001332192795501),
        ]
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.0, 0.0),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tbsv('u', 't', 'n', 6, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-7.6101460816946579E-003, -9.4660182789718411E-002),
            Complex::new(6.8147857987214427E-002, 1.5234806158371664E-003),
            Complex::new(2.9528998152503227, -0.17560477971933142),
            Complex::new(3.2401253786020834, -6.2781243245851357),
            Complex::new(-2.6978707277566238, 7.2418066105756456),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tbsv('u', 't', 'u', 6, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.64901006221771240, 0.77214217185974121),
            Complex::new(-0.90450126675144427, -0.75266294427094904),
            Complex::new(-0.17021002669464957, -0.33332994457346288),
            Complex::new(-1.3114892067913817, -1.0467456831262383),
            Complex::new(2.3201126210298955, 1.3473362542512941),
            Complex::new(0.56213947293563105, -0.67778084986020015),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tbsv('l', 't', 'n', 6, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.26783294222237042, 5.3412212596831878),
            Complex::new(-2.4708523091856045, -10.504946494218570),
            Complex::new(-2.5360982977762716, -1.6617901325948061),
            Complex::new(7.8824301124023846, -2.4955417126761907),
            Complex::new(2.7655458811932285, -0.55943960377267044),
            Complex::new(0.20166478159136628, -4.9273708864892143E-002)
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tbsv('l', 't', 'u', 6, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-2.8778446406669276, -3.9777722375573310),
            Complex::new(-3.1812735166230288, -1.3246826397301936),
            Complex::new(0.87222461178592359, 0.17438318677524328),
            Complex::new(0.34197712105112110, 1.3634675771268798),
            Complex::new(0.20611923584205027, 1.1345702378188081),
            Complex::new(-0.11775359511375427, -0.27577802538871765),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tbsv('u', 'c', 'n', 6, 3, &a, 6, &mut x, -1);
    capproximately!(
        x,
        vec![
            Complex::new(-2.4553590852142104, -5.3900057558999848),
            Complex::new(3.4605613317430035, 3.7946358385801675),
            Complex::new(2.3680417148052602, 6.5876435306548492E-002),
            Complex::new(0.36170019016395960, -8.2795869074406506E-003),
            Complex::new(0.15114802344991249, 0.41203114326506957),
            Complex::new(-0.13368054624873069, -0.18738554063649901),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tbsv('u', 'c', 'u', 6, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.64901006221771240, 0.77214217185974121),
            Complex::new(0.89701536056711850, 0.76156909872986789),
            Complex::new(0.71752083895656837, 0.41283609238991770),
            Complex::new(1.1504109107461260, 1.4022465442979499),
            Complex::new(-0.36504790906346596, -1.4582194308823537),
            Complex::new(0.95516886987466565, -0.74130583244324910),
        ]
    );
    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tbsv('l', 'c', 'n', 6, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(1.2134872328945268, 3.3689875258937509),
            Complex::new(-0.51454080049378648, -11.238018851608174),
            Complex::new(1.1080551530027434, -1.5280652728327344),
            Complex::new(-7.9545835289539690, -7.6831586104367161E-002),
            Complex::new(-2.6639255157882884, 0.22590102372385532),
            Complex::new(-0.17505966769138140, 0.11158268354623273)
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tbsv('l', 'c', 'u', 6, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-2.1629786169498084, -2.8234267947755658),
            Complex::new(2.6661358253363017, -3.4888111144575142),
            Complex::new(0.60070653051319522, 0.25805634910918834),
            Complex::new(0.10717711933379223, 1.3226800014067428),
            Complex::new(0.20002703546058132, 1.1371715273742764),
            Complex::new(-0.11775359511375427, -0.27577802538871765),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(0.0, 0.0),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(0.0, 0.0),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tbsv('l', 'c', 'u', 0, 3, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.64901006221771240, 0.77214217185974121),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.66413569450378418, -0.42481029033660889),
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.14377148449420929, 0.99698686599731445),
            Complex::new(-0.11775359511375427, -0.27577802538871765),
        ]
    );

    let result = std::panic::catch_unwind(|| {
        complex::tbsv('x', 't', 'u', 6, 5, &a, 6, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::tbsv('u', 'x', 'u', 6, 5, &a, 6, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::tbsv('u', 't', 'c', 6, 5, &a, 6, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::tbsv('u', 't', 'n', -6, 5, &a, 6, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::tbsv('u', 't', 'n', 6, -5, &a, 6, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::tbsv('u', 't', 'n', 6, 5, &a, 4, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::tbsv('u', 't', 'n', 6, 5, &a, 4, &mut vec![], 0);
    });
    assert!(result.is_err());
}

#[test]
fn tpmv() {
    let mut ap = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut ap, 6, 6, 0.0);
    let ap = matrix::complex::pack_upper(ap, 6, 6, 5);
    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tpmv('u', 'n', 'n', 6, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.21480810469340378, 0.69579227075292904),
            Complex::new(-0.26094969695022141, 2.1899200001689807),
            Complex::new(-0.33737654149548413, 0.95849359062438277),
            Complex::new(0.27866515796963282, -3.2357536489400593),
            Complex::new(-0.44083150014179662, -1.4766224545245905),
            Complex::new(-0.13564174811293128, -0.31767194384784148),
        ]
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tpmv('u', 'n', 'u', 6, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.21480810469340378, 0.69579227075292904),
            Complex::new(-0.15915947579584566, 1.7666898656354610),
            Complex::new(0.68550578102442095, 0.65675536979746663),
            Complex::new(1.8959032428106883, -1.0481669938732936),
            Complex::new(0.12916224198601700, 0.76626246552352040),
            Complex::new(-0.11775359511375427, -0.27577802538871765),
        ]
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tpmv('l', 'n', 'n', 6, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(2.9900855818740801E-003, -0.12544832535154660),
            Complex::new(0.21805658972259057, -0.64067279808502109),
            Complex::new(0.55454303196612764, 3.3278109254821664E-002),
            Complex::new(-1.9443462712775996, -0.20935680748047414),
            Complex::new(-1.4832770655398049, 1.4177745116425551),
        ]
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tpmv('l', 'n', 'u', 6, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.11916876584291458, -0.21951562166213989),
            Complex::new(0.94732568240280290, -0.22682655495516757),
            Complex::new(1.7616956826635151, -0.88734821132789865),
            Complex::new(-1.6881006148507482, 1.4323906112766398),
            Complex::new(-1.4653889125406279, 1.4596684301016789),
        ]
    );
    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    for c in ap.iter() {
        println!("{:?}", c.im);
    }
    complex::tpmv('u', 't', 'n', 6, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.22095898699729033, 0.20371451287137976),
            Complex::new(-0.26583054006009044, -9.1245991110092461E-002),
            Complex::new(-0.45535556765459351, -2.4279571343292985),
            Complex::new(-1.1230038346504205, -2.3902810640005123),
            Complex::new(-1.1355567240741828, 1.0540975185653942),
        ]
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tpmv('u', 't', 'u', 6, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.11916876584291458, -0.21951562166213989),
            Complex::new(0.75705178245981464, -0.39298421193700861),
            Complex::new(1.1618825171864620, -0.24037047926253274),
            Complex::new(-0.55301009252260691, -0.14739614395240164),
            Complex::new(-1.1176685710750058, 1.0959914370245181),
        ]
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tpmv('l', 't', 'n', 6, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-2.4976814036793784, -0.23946010250133121),
            Complex::new(1.3996314833238857, 0.78835332165500283),
            Complex::new(2.0129874112512680, -1.9013967150345032),
            Complex::new(-0.36287068806792244, 0.10791794560781565),
            Complex::new(-0.12708341444083437, -0.87548495323359354),
            Complex::new(-0.13564174811293128, -0.31767194384784148)
        ]
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tpmv('l', 't', 'u', 6, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-2.4976814036793784, -0.23946010250133121),
            Complex::new(1.2774726318990970, 0.69428602534440953),
            Complex::new(2.7422565039314803, -1.4875504719046497),
            Complex::new(0.84428196262946487, -0.81270837497490467),
            Complex::new(0.12916224198601700, 0.76626246552352040),
            Complex::new(-0.11775359511375427, -0.27577802538871765),
        ]
    );
    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tpmv('u', 'c', 'n', 6, &ap, &mut x, -1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.99053756426170603, -0.22489208505092062),
            Complex::new(0.37551077322001747, 1.4700155236647587),
            Complex::new(1.3389624021526658, 0.18727475284257045),
            Complex::new(0.24995649269782261, 0.50298068683844832),
            Complex::new(-0.86662654194144162, 3.1989469298911999E-002),
            Complex::new(-3.0267127927761095E-002, -0.39887173640357432),
        ]
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tpmv('u', 'c', 'u', 6, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.11916876584291458, -0.21951562166213989),
            Complex::new(0.64020916573379671, -0.32955369010140600),
            Complex::new(1.0371188810777623, 0.13627095379948040),
            Complex::new(0.64128318005638185, 1.5767700095488548),
            Complex::new(-0.78512600062635252, 0.95230004796124257),
        ]
    );
    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tpmv('l', 'c', 'n', 6, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.21231009546066959, 1.4470766140474218),
            Complex::new(-3.5191155812410297, 1.0437484124731138),
            Complex::new(1.1505912212692970, -0.10552584409884025),
            Complex::new(-0.21019657406961945, -0.75391528212598224),
            Complex::new(-0.23082261213115096, -0.82038820112026123),
            Complex::new(-0.13564174811293128, -0.31767194384784148),
        ]
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tpmv('l', 'c', 'u', 6, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.21231009546066959, 1.4470766140474218),
            Complex::new(-3.5314545767923575, 0.89006305075836589),
            Complex::new(1.0806654353856708, -0.94111695640190396),
            Complex::new(1.3036001790457676, -0.86876185148143970),
            Complex::new(-1.2774295761460941E-002, 0.82686750312289536),
            Complex::new(-0.11775359511375427, -0.27577802538871765),
        ]
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tpmv('l', 'c', 'u', 0, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.11916876584291458, -0.21951562166213989),
            Complex::new(0.66413569450378418, -0.42481029033660889),
            Complex::new(1.1009690761566162, -0.41898009181022644),
            Complex::new(0.14377148449420929, 0.99698686599731445),
            Complex::new(-0.11775359511375427, -0.27577802538871765)
        ]
    );

    let result = std::panic::catch_unwind(|| {
        complex::tpmv('x', 't', 'u', 6, &ap, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::tpmv('u', 'x', 'u', 6, &ap, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::tpmv('u', 't', 'x', 6, &ap, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::tpmv('u', 't', 'n', -6, &ap, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::tpmv('u', 't', 'n', 6, &ap, &mut vec![], 0);
    });
    assert!(result.is_err());
}

#[test]
fn tpsv() {
    let mut ap = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut ap, 6, 6, 0.0);
    let ap = matrix::complex::pack_upper(ap, 6, 6, 5);
    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(0.0, 0.0),
    ];
    complex::tpsv('u', 'n', 'n', 6, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-3.3851017924301425, -1.3820753807292496),
            Complex::new(-1.3324989852266680, 0.67573580100024855),
            Complex::new(-2.0061115606914002, 0.30873159135263567),
            Complex::new(0.53366394452817878, 0.69694725868851992),
            Complex::new(3.3264201117576032E-002, -0.76983395647848252),
            Complex::new(0.0000000000000000, 0.0000000000000000),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(0.0, 0.0),
    ];
    complex::tpsv('u', 'n', 'u', 6, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.15547987558715703, -2.1970472032566155),
            Complex::new(-0.64651555475774991, -0.84393783126494526),
            Complex::new(0.41271091050372538, -1.3369077805429912),
            Complex::new(0.38779544173476843, 0.38706108295362363),
            Complex::new(0.14377148449420929, 0.99698686599731445),
            Complex::new(0.0000000000000000, 0.0000000000000000),
        ]
    );
    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tpsv('l', 'n', 'n', 6, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.42327154210319390, -0.26082687970715562),
            Complex::new(0.39222678598002386, -0.23224612883036766),
            Complex::new(-1.2997678064087332, -0.94478598552509907),
            Complex::new(-0.61517805667934433, -0.20792193349925353),
            Complex::new(6.7048390491203161E-002, -0.61100462807805289),
        ]
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tpsv('l', 'n', 'u', 6, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.11916876584291458, -0.21951562166213989),
            Complex::new(0.38094570660476545, -0.62279402571805020),
            Complex::new(0.77065302766175314, 0.23155213977648836),
            Complex::new(1.9235986450552391, 0.62921553153130372),
            Complex::new(-0.67065782778448713, -2.0842393309515446),
        ]
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tpsv('u', 't', 'n', 6, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.20114145280773280, -5.1337030960433162E-002),
            Complex::new(-0.11347507554592454, 1.7911796556818167),
            Complex::new(0.17118012099267396, 0.41319272538503060),
            Complex::new(0.11383486554884616, 0.27127658857574360),
            Complex::new(0.13250679287302938, 0.40677081363532691),
        ]
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tpsv('u', 't', 'u', 6, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.11916876584291458, -0.21951562166213989),
            Complex::new(0.57121960654775372, -0.45663636873620916),
            Complex::new(1.0247114696420028, -0.62076677198172869),
            Complex::new(1.0568841164395277, 1.9635907758261411),
            Complex::new(0.37156836221836176, -2.7223246591641685),
        ]
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tpsv('l', 't', 'n', 6, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(12.179271994541113, 14.308616957116691),
            Complex::new(20.022860582238643, -11.394470245205239),
            Complex::new(5.3345587897852154, 2.9774530389843976),
            Complex::new(-3.4438036297245773, 0.90496225731641267),
            Complex::new(-0.18635839568165291, -1.8489702849007796),
            Complex::new(-0.10222449470844058, -0.23940899018683348),
        ],
        1E-5
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tpsv('l', 't', 'u', 6, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-4.6890622533536845, -4.3105166275435485),
            Complex::new(-4.9095388171859451, 1.5119249268050328),
            Complex::new(-2.0011605096316609, 0.28698531071705224),
            Complex::new(1.4229022368946207, 9.6376644440676285E-002),
            Complex::new(0.15838072700240158, 1.2277112664711085),
            Complex::new(-0.11775359511375427, -0.27577802538871765),
        ]
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tpsv('u', 'c', 'n', 6, &ap, &mut x, -1);
    capproximately!(
        x,
        vec![
            Complex::new(0.39803859497384503, -1.0008794896314415),
            Complex::new(-0.95461016818598643, -0.10455630604512522),
            Complex::new(0.10819566776526710, -0.48235598200894436),
            Complex::new(-2.1147145699491019, -1.3663776876599896),
            Complex::new(0.80628269206021930, -8.8111846075539443E-002),
            Complex::new(-0.15013342830572224, -0.16730119413206657),
        ]
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tpsv('u', 'c', 'u', 6, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.11916876584291458, -0.21951562166213989),
            Complex::new(0.68806222327377164, -0.52006689057181177),
            Complex::new(1.1537416954878212, -0.99972441484002239),
            Complex::new(-0.82757703535350080, 7.3583568030202828E-002),
            Complex::new(1.5551601886103397, -1.4986169107249498),
        ]
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tpsv('l', 'c', 'n', 6, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-16.223221775482749, 12.881420892412462),
            Complex::new(-14.487432918075610, -15.055884393595845),
            Complex::new(-6.4169359859497845, -0.96002124817834023),
            Complex::new(1.5382948185317014, 4.3925279091258904),
            Complex::new(-0.48218019458408296, -1.7482432429647807),
            Complex::new(-0.10222449470844058, -0.23940899018683348),
        ],
        1E-5
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tpsv('l', 'c', 'u', 6, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(2.5181705640049294, -5.4011920508919671),
            Complex::new(3.8304845810211443, -1.6707772413094739),
            Complex::new(0.73533020104709457, -0.53087771067682288),
            Complex::new(0.94106504084051901, 0.16204542951614687),
            Complex::new(0.30031726474987952, 1.1671062288717335),
            Complex::new(-0.11775359511375427, -0.27577802538871765),
        ]
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::tpsv('l', 'c', 'u', 0, &ap, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.11916876584291458, -0.21951562166213989),
            Complex::new(0.66413569450378418, -0.42481029033660889),
            Complex::new(1.1009690761566162, -0.41898009181022644),
            Complex::new(0.14377148449420929, 0.99698686599731445),
            Complex::new(-0.11775359511375427, -0.27577802538871765),
        ]
    );

    let result = std::panic::catch_unwind(|| {
        complex::tpsv('x', 't', 'u', 6, &ap, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::tpsv('u', 'x', 'u', 6, &ap, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::tpsv('u', 't', 'x', 6, &ap, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::tpsv('u', 't', 'n', -6, &ap, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::tpsv('u', 't', 'n', 6, &ap, &mut vec![], 0);
    });
    assert!(result.is_err());
}

#[test]
fn trmv() {
    let mut a = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 6), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut a, 6, 6, 0.0);
    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(0.0, 0.0),
    ];
    complex::trmv('u', 'n', 'n', 6, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.83150675414885056, 1.9278190152391241),
            Complex::new(-0.39338349525210226, 2.0694548243325235),
            Complex::new(-0.40640543562459586, 0.84668404121139496),
            Complex::new(0.19690462573740852, -3.4126079216408423),
            Complex::new(-0.42622225763360433, -1.2458980540507965),
            Complex::new(0.0000000000000000, 0.0000000000000000),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(0.0, 0.0),
    ];
    complex::trmv('u', 'n', 'u', 6, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.99249198509378145, 1.4460226393835436),
            Complex::new(-0.29159327409772651, 1.6462246897990038),
            Complex::new(0.61647688689530922, 0.54494582038447881),
            Complex::new(1.8141427105784640, -1.2250212665740765),
            Complex::new(0.14377148449420929, 0.99698686599731445),
            Complex::new(0.0000000000000000, 0.0000000000000000),
        ]
    );

    let mut a = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 6), 6, 1, 6, 1, 6);
    matrix::complex::set_upper(&mut a, 6, 6, 0.0);

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(0.0, 0.0),
    ];
    complex::trmv('l', 'n', 'n', 6, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.48802483127278151, 1.2539385477153218),
            Complex::new(-0.96537710222919570, -0.85185586585022488),
            Complex::new(-1.2393312224660826, 1.2132717007075515),
            Complex::new(-2.5654559956894172, -1.3102608123698012),
            Complex::new(-1.5184106544536302, -3.6775354290217734),
            Complex::new(3.4744621149400827, -0.57679995107873072)
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(0.0, 0.0),
    ];
    complex::trmv('l', 'n', 'u', 6, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.64901006221771240, 0.77214217185974121),
            Complex::new(-0.86358688107481996, -1.2750860003837445),
            Complex::new(-0.21644889994617750, 0.91153347988063538),
            Complex::new(-0.94821791084836171, 0.87732584269696456),
            Complex::new(-0.94841691232581660, -1.4346505089736628),
            Complex::new(3.4744621149400827, -0.57679995107873072),
        ]
    );

    let mut a = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 6), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut a, 6, 6, 0.0);
;
    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(0.0, 0.0),
    ];
    complex::trmv('u', 't', 'n', 6, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.48802483127278151, 1.2539385477153218),
            Complex::new(1.0241460862502108, 2.6733200210011354E-002),
            Complex::new(0.65454049159191263, -0.82986178215861273),
            Complex::new(-0.93126308391969648, -2.2538919124303631),
            Complex::new(-1.0774122839169915, -2.4272069861881223),
            Complex::new(-1.1294976357766762, 0.95227101262259106),
        ]
    );

    let mut a = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 6), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut a, 6, 6, 0.0);
;
    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(0.0, 0.0),
    ];
    complex::trmv('u', 't', 'u', 6, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.64901006221771240, 0.77214217185974121),
            Complex::new(1.1259363074045865, -0.39649693432350830),
            Complex::new(1.6774228141118177, -1.1316000029855289),
            Complex::new(0.68597500092135899, -6.6305257363597381E-002),
            Complex::new(-0.50741854178917789, -0.18432206614001156),
            Complex::new(-1.1294976357766762, 0.95227101262259106),
        ]
    );

    let mut a = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 6), 6, 1, 6, 1, 6);
    matrix::complex::set_upper(&mut a, 6, 6, 0.0);

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(0.0, 0.0),
    ];
    complex::trmv('l', 't', 'n', 6, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(2.2251660077830033, 2.3238341766138930),
            Complex::new(0.27152304192522392, -2.3031402183418956),
            Complex::new(1.0402809982832415, 2.8743260159852113),
            Complex::new(0.76270793503965795, -2.6553038624241152),
            Complex::new(-0.42622225763360433, -1.2458980540507965),
            Complex::new(0.0000000000000000, 0.0000000000000000),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(0.0, 0.0),
    ];
    complex::trmv('l', 't', 'u', 6, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(2.0641807768380724, 1.8420378007583125),
            Complex::new(0.37331326307959944, -2.7263703528754153),
            Complex::new(2.0631633208031466, 2.5725877951582952),
            Complex::new(2.3799460198807134, -0.46771720735734945),
            Complex::new(0.14377148449420929, 0.99698686599731445),
            Complex::new(0.0000000000000000, 0.0000000000000000),
        ]
    );

    let mut a = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 6), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut a, 6, 6, 0.0);
;
    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(0.0, 0.0),
    ];
    complex::trmv('u', 'c', 'n', 6, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-1.1513151820979886, 0.69642190434815632),
            Complex::new(0.25139557645729682, -1.3313132788431306),
            Complex::new(0.50669212996350743, -0.56114501912958303),
            Complex::new(1.1938025271467749, 3.3443289488025765),
            Complex::new(0.58289064685587677, -0.78703924639420531),
            Complex::new(-0.23183789661954024, 1.2835785838419986),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(0.0, 0.0),
    ];
    complex::trmv('u', 'c', 'u', 6, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.64901006221771240, 0.77214217185974121),
            Complex::new(-0.15897511885275151, -1.4765058693244129),
            Complex::new(1.2095200637738619, -1.3632466929137586),
            Complex::new(0.94750069304083251, 0.63502464640267542),
            Complex::new(0.66981737160638177, 1.5255068342641325),
            Complex::new(-0.23183789661954024, 1.2835785838419986),
        ]
    );

    let mut a = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 6), 6, 1, 6, 1, 6);
    matrix::complex::set_upper(&mut a, 6, 6, 0.0);

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(0.0, 0.0),
    ];
    complex::trmv('l', 'c', 'n', 6, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.90061994398728085, -1.5995329267787675),
            Complex::new(5.3055233216601705, 1.9450131306691401),
            Complex::new(-2.2713161685357997, -1.7723500578825266),
            Complex::new(0.10663390086878710, 2.6049305359109640),
            Complex::new(5.6844759743704287E-002, -1.3155592146610233),
            Complex::new(0.0000000000000000, 0.0000000000000000),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(0.0, 0.0),
    ];
    complex::trmv('l', 'c', 'u', 6, &a, 6, &mut x, -1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.64901006221771240, 0.77214217185974121),
            Complex::new(-1.2498731094504563, -1.1083395882641662),
            Complex::new(0.68069188988759555, 0.24910421006826589),
            Complex::new(0.47157803291335676, -2.7400158851063168),
            Complex::new(3.6613593566644753, 2.0828881554261853),
            Complex::new(3.3443933864997035, -3.9680368366347851),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(0.0, 0.0),
    ];
    complex::trmv('l', 'c', 'u', 0, &a, 6, &mut x, -1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.64901006221771240, 0.77214217185974121),
            Complex::new(-0.11916876584291458, -0.21951562166213989),
            Complex::new(0.66413569450378418, -0.42481029033660889),
            Complex::new(1.1009690761566162, -0.41898009181022644),
            Complex::new(0.14377148449420929, 0.99698686599731445),
            Complex::new(0.0000000000000000, 0.0000000000000000),
        ]
    );

    let result = std::panic::catch_unwind(|| {
        complex::trmv('x', 't', 'u', 6, &a, 6, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::trmv('u', 'x', 'u', 6, &a, 6, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::trmv('u', 't', 'x', 6, &a, 6, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::trmv('u', 't', 'n', -6, &a, 6, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::trmv('u', 't', 'n', 6, &a, 6, &mut vec![], 0);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::trmv('u', 't', 'n', 6, &a, 5, &mut vec![], 1);
    });
    assert!(result.is_err());
}

#[test]
fn trsv() {
    let mut a = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut a, 6, 6, 0.0);
;
    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(0.0, 0.0),
    ];
    complex::trsv('u', 'n', 'n', 6, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-3.3851017924301425, -1.3820753807292496),
            Complex::new(-1.3324989852266680, 0.67573580100024855),
            Complex::new(-2.0061115606914002, 0.30873159135263567),
            Complex::new(0.53366394452817878, 0.69694725868851992),
            Complex::new(3.3264201117576032E-002, -0.76983395647848252),
            Complex::new(0.0000000000000000, 0.0000000000000000),
        ]
    );
    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(0.0, 0.0),
    ];
    complex::trsv('u', 'n', 'u', 6, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.15547987558715703, -2.1970472032566155),
            Complex::new(-0.64651555475774991, -0.84393783126494526),
            Complex::new(0.41271091050372538, -1.3369077805429912),
            Complex::new(0.38779544173476843, 0.38706108295362363),
            Complex::new(0.14377148449420929, 0.99698686599731445),
            Complex::new(0.0000000000000000, 0.0000000000000000),
        ]
    );

    let mut a = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_upper(&mut a, 6, 6, 0.0);

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::trsv('l', 'n', 'n', 6, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.20114145280773280, -5.1337030960433162E-002),
            Complex::new(-0.58779672483646239, 1.4444342253366556),
            Complex::new(-0.59996084429980345, 1.8923968370248347),
            Complex::new(2.3943121537219199, 1.2601259090017753),
            Complex::new(-1.1716441125457138, 2.4468331626047548),
        ]
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::trsv('l', 'n', 'u', 6, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.11916876584291458, -0.21951562166213989),
            Complex::new(0.89736221505014213, -0.55306139391414555),
            Complex::new(0.74985303977344109, -2.4958541378702988),
            Complex::new(3.2482401623515491, 3.4359461997861374),
            Complex::new(-6.1115892109325278, 5.4105978948754974),
        ]
    );

    let mut a = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut a, 6, 6, 0.0);
;
    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::trsv('u', 't', 'n', 6, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(0.20114145280773280, -5.1337030960433162E-002),
            Complex::new(-0.11347507554592454, 1.7911796556818167),
            Complex::new(0.17118012099267396, 0.41319272538503060),
            Complex::new(0.11383486554884616, 0.27127658857574360),
            Complex::new(0.13250679287302938, 0.40677081363532691),
        ]
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::trsv('u', 't', 'u', 6, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.11916876584291458, -0.21951562166213989),
            Complex::new(0.57121960654775372, -0.45663636873620916),
            Complex::new(1.0247114696420028, -0.62076677198172869),
            Complex::new(1.0568841164395277, 1.9635907758261411),
            Complex::new(0.37156836221836176, -2.7223246591641685),
        ]
    );

    let mut a = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_upper(&mut a, 6, 6, 0.0);

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(0.0, 0.0),
    ];
    complex::trsv('l', 't', 'n', 6, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(2.5983253692194772, -15.102300459040050),
            Complex::new(6.3141701681529163, -6.9203837691759773),
            Complex::new(-3.7488703711821940, 6.4167258382816481),
            Complex::new(0.27337260555960724, 0.88409421189284909),
            Complex::new(3.3264201117576032E-002, -0.76983395647848252),
            Complex::new(0.0000000000000000, 0.0000000000000000),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(0.0, 0.0),
    ];
    complex::trsv('l', 't', 'u', 6, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.42363350716824644, -2.2093867649350898),
            Complex::new(2.2225892101670914, -1.7538965376185174),
            Complex::new(-1.1422197242316745, -0.27970358425951164),
            Complex::new(-0.17800786756748099, -0.37024297626310343),
            Complex::new(0.14377148449420929, 0.99698686599731445),
            Complex::new(0.0000000000000000, 0.0000000000000000),
        ]
    );
    let mut a = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_lower(&mut a, 6, 6, 0.0);
;
    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::trsv('u', 'c', 'n', 6, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.15262265513448822, 0.14071140652106690),
            Complex::new(-1.4443838535530182, -0.39906073251810770),
            Complex::new(-4.6332350725494867E-002, -0.20224298804899468),
            Complex::new(-1.4838666079245511, -0.94099706199515598),
            Complex::new(0.39843405102446616, -0.29059808376945545),
        ]
    );

    let mut x = vec![
        Complex::new(0.0, 0.0),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(-0.11775359816595128, -0.27577802908802723),
    ];
    complex::trsv('u', 'c', 'u', 6, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(0.0000000000000000, 0.0000000000000000),
            Complex::new(-0.11916876584291458, -0.21951562166213989),
            Complex::new(0.68806222327377164, -0.52006689057181177),
            Complex::new(1.1537416954878212, -0.99972441484002239),
            Complex::new(-0.82757703535350080, 7.3583568030202828E-002),
            Complex::new(1.5551601886103397, -1.4986169107249498),
        ]
    );

    let mut a = matrix::complex::slice(fixtures::complex::matrix_mxn(6, 8), 6, 1, 6, 1, 6);
    matrix::complex::set_upper(&mut a, 6, 6, 0.0);

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(0.0, 0.0),
    ];
    complex::trsv('l', 'c', 'n', 6, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(1.4997211511167985, 1.6317882003200239),
            Complex::new(0.66507043863638171, -0.28918592454053232),
            Complex::new(-0.57313915841867680, -0.55599486430937639),
            Complex::new(1.2909233609283014E-002, -9.3569865927350199E-002),
            Complex::new(-0.24941512573253091, -0.72906990246418102),
            Complex::new(0.0000000000000000, 0.0000000000000000),
        ]
    );

    let mut x = vec![
        Complex::new(-0.6490100777088978, 0.7721421858045301),
        Complex::new(-0.11916876241803812, -0.21951562675343952),
        Complex::new(0.6641356998941105, -0.4248102833772871),
        Complex::new(1.100969102194087, -0.418980099421959),
        Complex::new(0.14377148075806995, 0.9969868609091059),
        Complex::new(0.0, 0.0),
    ];
    complex::trsv('l', 'c', 'u', 6, &a, 6, &mut x, -1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.64901006221771240, 0.77214217185974121),
            Complex::new(1.0115355777646271, 0.66930834493988645),
            Complex::new(1.6199462693950109, -2.6460116567259719),
            Complex::new(7.9681750368777013, 2.1980136848408875),
            Complex::new(-8.2087922632419357, -2.9159109438871607),
            Complex::new(-7.9036649818434919, -8.2778612243653438),
        ]
    );

    let mut x = vec![
        Complex::new(-0.64901006221771240, 0.77214217185974121),
        Complex::new(-0.11916876584291458, -0.21951562166213989),
        Complex::new(0.66413569450378418, -0.42481029033660889),
        Complex::new(1.1009690761566162, -0.41898009181022644),
        Complex::new(0.14377148449420929, 0.99698686599731445),
        Complex::new(-0.11775359511375427, -0.27577802538871765),
    ];
    complex::trsv('l', 'c', 'u', 0, &a, 6, &mut x, 1);
    capproximately!(
        x,
        vec![
            Complex::new(-0.64901006221771240, 0.77214217185974121),
            Complex::new(-0.11916876584291458, -0.21951562166213989),
            Complex::new(0.66413569450378418, -0.42481029033660889),
            Complex::new(1.1009690761566162, -0.41898009181022644),
            Complex::new(0.14377148449420929, 0.99698686599731445),
            Complex::new(-0.11775359511375427, -0.27577802538871765),
        ]
    );

    let result = std::panic::catch_unwind(|| {
        complex::trsv('x', 't', 'u', 6, &a, 6, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::trsv('u', 'x', 'u', 6, &a, 6, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::trsv('u', 't', 'x', 6, &a, 6, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::trsv('u', 't', 'n', -6, &a, 6, &mut vec![], 1);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::trsv('u', 't', 'n', 6, &a, 6, &mut vec![], 0);
    });
    assert!(result.is_err());

    let result = std::panic::catch_unwind(|| {
        complex::trsv('u', 't', 'n', 6, &a, 5, &mut vec![], 1);
    });
    assert!(result.is_err());
}
