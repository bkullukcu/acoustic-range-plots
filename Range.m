% -------------------------------------------------------------------
%  Generated by MATLAB on 6-Nov-2020 21:51:59
%  MATLAB version: 9.8.0.1359463 (R2020a) Update 1
% -------------------------------------------------------------------
saveVarsMat = load('Range.mat');

B = [2.7604733097482447E-5 0.0020392017774081977 0.0076748852606828028 ...
     0.017065212042668025 0.030210182560940494 0.047109796869469873 0.067764054981579142 ...
     0.092172956901834455 0.12033650263214343 0.15225469217341761 0.18792752552613695 ...
     0.22735500269057382 0.27053712366689225 0.31747388845519547 0.36816529705555134 ...
     0.42261134946800616 0.48081204569259162 0.54276738572933114 0.60847736957824106 ...
     0.677941997239334 0.75116126871261923 0.82813518399810437 0.90886374309579465 ...
     0.99334694600569484 1.0815847927278079 1.1735772832621378 1.2693244176086858 ...
     1.3688261957674541 1.4720826177384438 1.5790936835216578 1.6898593931170949 ...
     1.8043797465247577 1.922654743744646 2.0446843847767617 2.1704686696211044 ...
     2.300007598277674 2.4333011707464722 2.5703493870274987 2.7111522471207543 ...
     2.8557097510262377 3.0040218987439515 3.1560886902738932 3.3119101256160661 ...
     3.4714862047704687 3.6348169277371003 3.8019022945159624 3.9727423051070545 ...
     4.147336959510378 4.3256862577259305 4.507790199753714 4.6936487855937274 ...
     4.8832620152459718 5.0766298887104462 5.2737524059871523 5.4746295670760912 ...
     5.6792613719772573 5.8876478206906553 6.0997889132162859 6.3156846495541448 ...
     6.5353350297042345 6.7587400536665561 6.985899721441112 7.2168140330278963 ...
     7.4514829884269105 7.6899065876381556 7.9320848306616352 8.1780177174973439 ...
     8.4277052481452817 8.6811474226054557 8.9383442408778553 9.1992957029624876 ...
     9.46400180885935 9.7324625585684483 10.004677952089775 10.280647989423333 ...
     10.560372670569119 10.843851995527144 11.131085964297394 11.422074576879877 ...
     11.716817833274593 12.015315733481538 12.317568277500715 12.623575465332122 ...
     12.933337296975765 13.246853772431631 13.564124891699732 13.885150654780066 ...
     14.209931061672636 14.538466112377426 14.870755806894454 15.206800145223715 ...
     15.546599127365205 15.890152753318926 16.237461023084872 16.588523936663062 ...
     16.943341494053477 17.301913695256125 17.664240540271 18.030322029098112 ...
     18.400158161737451];

B1 = [2.7474779987156006E-7 0.023500160187483032 0.031213012006026719 0.041035304317185259 ...
      0.054335437864507244 0.07130759307052463 0.0920014511049012 0.11643429762648279 ...
      0.14461341157505064 0.17654228712303321 0.21222277038704554 0.25165591100942558 ...
      0.29484234145102739 0.34178246118369537 0.39247653256052889 0.44692473363899127 ...
      0.50512718870579232 0.56708398664497273 0.63279519238259074 0.7022608542370854 ...
      0.77548100876742787 0.85245568404714778 0.93318490192231018 1.0176686795984293 ...
      1.1059070307749157 1.1978999664686987 1.2936474956206692 1.3931496255480056 ...
      1.4964063622855646 1.6034177108463852 1.7141836754224693 1.828704259541015 ...
      1.9469794661870214 2.0690092979003163 2.1947937568529041 2.3243328449110869 ...
      2.4576265636856771 2.5946749145728427 2.7354778987875403 2.8800355173910188 ...
      3.0283477713135816 3.1804146613734967 3.3362361882928222 3.4958123527106491 ...
      3.6591431551942986 3.8262285962487721 3.9970686763247905 4.1716633958256537 ...
      4.3500127551130934 4.5321167545123213 4.7179753943163565 4.9075886747897819 ...
      5.1009565961719883 5.29807915868 5.4989563625109419 5.7035882078441755 ...
      5.9119746948432006 6.1241158236573039 6.3400115944230109 6.55966200726538 ...
      6.7830670622991365 7.0102267596296874 7.2411410993539924 7.475810081561395 ...
      7.7142337063343049 7.9564119737488461 8.2023448838754174 8.4520324367792039 ...
      8.7054746325206445 8.962671471155808 9.2236229527368074 9.4883290773121054 ...
      9.7567898449268267 10.029005255623012 10.304975309439904 10.584700006414124 ...
      10.868179346579923 11.155413329969321 11.446401956612313 11.741145226537004 ...
      12.039643139769748 12.341895696335291 12.647902896256873 12.95766473955635 ...
      13.271181226254274 13.588452356370004 13.90947812992178 14.2342585469268 ...
      14.562793607401275 14.895083311360533 15.231127658819046 15.570926649790493 ...
      15.914480284287816 16.261788562323261 16.612851483908429 16.967669049054326 ...
      17.326241257771354 17.688568110069426 18.054649605957923 18.424485745445754 ...
      ];

B10 = [2.058998456912949E-8 0.015427761084878514 0.046457371311503529 0.089292223574243934 ...
       0.14120026397092941 0.19861523212773236 0.25870689449703338 0.31971713109938721 ...
       0.38080240932110276 0.44174693572994539 0.50270346880887 0.56400776693551524 ...
       0.6260627456738026 0.689273619545252 0.75401571722779648 0.82062144714276242 ...
       0.8893776469402157 0.96052806894944609 1.0342780433648149 1.1107997454665062 ...
       1.1902372890788977 1.2727113069522615 1.358322909634766 1.4471570287248032 ...
       1.5392852013199187 1.634767870178276 1.7336562752961364 1.8359940063115852 ...
       1.941818276098406 2.0511609664751131 2.1640494882076426 2.2805074898467379 ...
       2.4005554434937193 2.5242111302562691 2.6514900438019171 2.7824057268897739 ...
       2.9169700529170508 3.0551934622294343 3.1970851611064348 3.3426532898557224 ...
       3.4919050652624035 3.6448469016820062 3.8014845142934677 3.9618230074033671 ...
       4.1258669501858787 4.2936204418307753 4.4650871677357973 4.6402704481049843 ...
       4.8191732800892737 5.00179837442042 5.1881481873365249 5.3782249484711109 ...
       5.5720306852728481 5.7695672444359447 5.9708363107484459 6.1758394237049377 ...
       6.3845779921791879 6.5970533074093565 6.8132665545122917 7.0332188227129961 ...
       7.2569111144493386 7.484344353490247 7.71551939218692 7.9504370179606187 ...
       8.18909795911697 8.43150289006509 8.67765243600975 8.9275471771762 ...
       9.181187652619851 9.4385743636665183 9.6997077770234483 9.9645883275964664 ...
       10.233216421044411 10.505592436098333 10.781716726669869 11.06158962377021 ...
       11.345211437258895 11.632582457439289 11.923702956515974 12.2185731899274 ...
       12.517193397565821 12.819563804895321 13.125684623977397 13.435556054412796 ...
       13.749178284207241 14.066551490568042 14.387675840637732 14.712551492170391 ...
       15.041178594155671 15.373557287395126 15.709687705034897 16.04956997305851 ...
       16.393204210743175 16.740590531082567 17.091729041179 17.446619842607276 ...
       17.805263031752784 18.167658700125671 18.53380693465315 18.903707817951577 ...
       ];

B2 = [1.1548050853248958E-7 0.035295356073316464 0.0609405734165841 0.077228020809831072 ...
      0.093248531167380244 0.11157800492546115 0.13303754687745348 0.15794217465691718 ...
      0.18643166828090096 0.21857527145334368 0.25441032403739389 0.29395835693488992 ...
      0.33723247527890915 0.38424101884046624 0.43498949626192573 0.48948166301735407 ...
      0.5477201498567642 0.60970684361528837 0.6754431258220418 0.7449300266518285 ...
      0.81816832685517116 0.89515862681255776 0.97590139428708378 1.0603969980623988 ...
      1.1486457320373566 1.2406478327490611 1.3364034922946235 1.4359128679817477 ...
      1.5391760896210087 1.6461932650959139 1.7569644846601835 1.8714898242840021 ...
      1.9897693482822889 2.1118031113957594 2.2375911604511605 2.3671335356951824 ...
      2.5004302718733156 2.63748139910789 2.7782869436169046 2.922846928305785 ...
      3.0711613732571186 3.2232302961379595 3.3790537125402134 3.5386316362663361 ...
      3.7019640795702071 3.8690510533609968 4.0398925673764037 4.2144886303304157 ...
      4.3928392500397617 4.5749444335325311 4.7608041871417424 4.9504185165862253 ...
      5.143787427040718 5.3409109231968293 5.541789009316159 5.7464216892767261 ...
      5.9548089666136805 6.16695084455503 6.3828473260531071 6.60249841381236 ...
      6.8259041103138953 7.0530644178372635 7.2839793384797717 7.5186488741737127 ...
      7.7570730267016863 7.99925179771031 8.24518518872247 8.4948732011483212 ...
      8.74831583629515 9.0055130953762443 9.2664649795188954 9.5311714897716069 ...
      9.7996326271106131 10.071848392445766 10.347818786625892 10.627543810443632 ...
      10.911023464639841 11.198257749907599 11.489246666895863 11.783990216212798 ...
      12.082488398428806 12.384741214079307 12.690748663667293 13.000510747665649 ...
      13.314027466519276 13.631298820647086 13.95232481044377 14.27710543628149 ...
      14.605640698511371 14.937930597464941 15.273975133455421 15.613774306778906 ...
      15.957328117715509 16.304636566530352 16.655699653474553 17.010517378786069 ...
      17.369089742690541 17.73141674540204 18.097498387123778 18.467334668048743 ...
      ];

B3 = [7.219684815051962E-8 0.0314578376549426 0.072995913494666886 0.10326504863595806 ...
      0.12747448871010039 0.1505608421668454 0.17493006584651302 0.20171531446006832 ...
      0.2314803293176301 0.26452470634189085 0.30101765566169764 0.34105993040869148 ...
      0.38471429451842648 0.43202138658616823 0.48300840778681708 0.53769409590501738 ...
      0.59609168619523378 0.65821073529377683 0.72405827927575328 0.79363958895848374 ...
      0.86695867447417274 0.94401862967099359 1.0248218717881632 1.109370311203481 ...
      1.1976654735896011 1.2897085891123936 1.3855006584391374 1.4850425021892009 ...
      1.5883347984026071 1.6953781112289852 1.8061729131088571 1.92071960207934 ...
      2.0390185153903588 2.1610699403028186 2.2868741227154552 2.4164312741049652 ...
      2.5497415771457344 2.6868051902884176 2.8276222515119751 2.97219288141526 ...
      3.1205171857776737 3.2725952576905151 3.4284271793393306 3.5880130235010035 ...
      3.7513528548065969 3.9184467308108824 4.08929470290163 4.2638968170755165 ...
      4.4422531146025381 4.6243636325968529 4.8102284045088259 4.999847460550451 ...
      5.1932208280642467 5.3903485318440865 5.5912305944149328 5.7958670362774232 ...
      6.0042578761222529 6.216403131018529 6.4323028165796554 6.6519569471097668 ...
      6.87536553573326 7.1025285945096206 7.3334461345354178 7.5681181660351022 ...
      7.8065446984419555 8.04872574047042 8.29466130018081 8.5443513850373716 ...
      8.797796001960382 9.0549951573730176 9.3159488572435958 9.5806571071236846 ...
      9.849119912182525 10.121337277238187 10.397309206785822 10.677035705023274 ...
      10.96051677587438 11.247752423010118 11.53874264986794 11.833487459669364 ...
      12.131986855436034 12.434240840004433 12.740249416039317 13.050012586046057 ...
      13.363530352381909 13.680802717266443 14.001829682791044 14.32661125092771 ...
      14.655147423537121 14.987438202376124 15.323483589104578 15.663283585291717 ...
      16.00683819242203 16.354147411900691 16.705211245058592 17.060029693156991 ...
      17.418602757391849 17.780930438897844 18.147012738752114 18.516849657977691 ...
      ];

B4 = [5.2532763666984274E-8 0.026579973905785471 0.071985701128335008 0.1142341572233098 ...
      0.149365905135242 0.18040216034119061 0.21030034731142622 0.24095313407762647 ...
      0.27348063227424318 0.30854653426606787 0.34655398105388913 0.38775558903430379 ...
      0.43231468798441325 0.48034002613267857 0.53190597529611217 0.58706464052438911 ...
      0.64585332700235865 0.70829926909850771 0.7744227005642712 0.844238893637504 ...
      0.91775954182964348 0.99499371567530881 1.075948534945194 1.1606296490592642 ...
      1.2490415855083554 1.3411880059915373 1.4370718970793084 1.536695713789578 ...
      1.6400614888709715 1.7471709168175327 1.8580254190601286 1.9726261949924682 ...
      2.0909742622350302 2.2130704886493406 2.3389156179752892 2.4685102905001322 ...
      2.6018550598278409 2.738950406566087 2.8797967495607337 3.0243944551667004 ...
      3.17274384493733 3.3248452020327708 3.4806987765852924 3.6403047902108705 ...
      3.8036634398186839 3.970774900840492 4.14163932997862 4.316256867552748 ...
      4.494627639511 4.6767517591590959 4.8626293286517521 5.0522604402829927 ...
      5.2456451776057031 5.4427836164057783 5.6436758255520187 5.8483218677395516 ...
      6.056721800141772 6.2688756749833754 6.4847835400452647 6.7044454391103976 ...
      6.9278614123583422 7.1550314967151509 7.38595572616428 7.6206341320233921 ...
      7.859066743191284 8.1012535863685287 8.3471946862550155 8.5968900657270915 ...
      8.8503397459967 9.1075437467545051 9.3685020862989337 9.6332147816525833 ...
      9.9016818486674563 10.1739033021202 10.449879155798437 10.729609422579125 ...
      11.013094114499768 11.300333242823182 11.591326818096588 11.886074850205432 ...
      12.184577348422581 12.486834321453307 12.79284577747646 13.1026117241822 ...
      13.416132168806589 13.733407118163418 14.054436578673407 14.379220556391092 ...
      14.707759057029593 15.040052085983479 15.376099648349816 15.715901748947635 ...
      16.059458392335944 16.406769582830364 16.757835324518574 17.112655621274584 ...
      17.471230476771989 17.833559894496315 18.199643877756415 18.569482429695121 ...
      ];

B5 = [4.1385674655231666E-8 0.023072239257954309 0.066735624437101848 0.11495449127130426 ...
      0.15930921144402377 0.19909875414506881 0.23613719957952004 0.27232522006833243 ...
      0.30909114088854761 0.34741426242797036 0.38795047110345371 0.43114010953645793 ...
      0.47728265961542432 0.52658529374524465 0.57919395022244258 0.63521329100136259 ...
      0.69471966981286448 0.75776968770830255 0.824405931233369 0.894660885442168 ...
      0.96855964576353559 1.0461218266712942 1.1273629248027075 1.2122953059097203 ...
      1.3009289286976504 1.3932718821195735 1.4893307887138536 1.5891111105856319 ...
      1.6926173838294369 1.7998533997932775 1.910822346457141 2.0255269196024588 ...
      2.1439694108982295 2.2661517782010767 2.3920757020425945 2.5217426313097056 ...
      2.6551538204100433 2.7923103596834276 2.9332132004224087 3.077863175563917 ...
      3.2262610168850094 3.3784073693600782 3.5343028032013986 3.6939478239995855 ...
      3.8573428812983592 4.0244883758733314 4.19538466593349 4.3700320724235011 ...
      4.5484308835724914 4.7305813588090944 4.9164837321414758 5.1061382150842061 ...
      5.299544999199953 5.4967042583128212 5.6976161504408216 5.9022808194874505 ...
      6.1106983967260575 6.3228690021054339 6.5387927454008237 6.7584697272308905 ...
      6.9819000399581528 7.2090837684879041 7.4400209909784474 7.6747117794737676 ...
      7.913156200468082 8.1553543154105927 8.4013061811574943 8.6510118503775129 ...
      8.9044713719162925 9.1616847911243351 9.4226521501526435 9.6873734882196043 ...
      9.9558488418522675 10.228078245104797 10.504061729756556 10.783799325491895 ...
      11.067291060063628 11.35453695944179 11.645537047949247 11.940291348385403 ...
      12.238799882139206 12.541062669292533 12.847079728714798 13.156851078149726 ...
      13.470376734294899 13.787656712874906 14.108691028708551 14.433479695770696 ...
      14.762022727249276 15.094320135597846 15.430371932584096 15.770178129334656 ...
      16.113738736376529 16.461053763675437 16.81212322067136 17.166947116311416 ...
      17.525525459080431 17.887858257029297 18.253945517801327 18.623787248656726 ...
      ];

B6 = [3.423002498804149E-8 0.020628401520738106 0.061147626618176167 0.11086084554279452 ...
      0.16101166143389545 0.20821363426268488 0.25246597010669974 0.29494113908495584 ...
      0.33693112649151719 0.37951292919268853 0.42350341112285361 0.4695005985291541 ...
      0.51793782749358463 0.56912959923481832 0.62330586628796725 0.68063648359464179 ...
      0.741248348569036 0.805237378518026 0.87267691368317368 0.94362365791474756 ...
      1.0181219187685 1.096206665461094 1.1779057579378271 1.2632415890741433 ...
      1.3522323071229625 1.4448927348434211 1.5412350672227535 1.6412694059851649 ...
      1.7450041726383743 1.8524464303009345 1.9636021364233569 2.0784763427182189 ...
      2.1970733544416037 2.3193968581382935 2.4454500247438844 2.5752355932981947 ...
      2.7087559393043272 2.8460131308526764 2.9870089749377082 3.1317450558691613 ...
      3.2802227672762796 3.4324433388928566 3.5884078590698478 3.7481172937740781 ...
      3.911572502683982 4.0787742528767392 4.2497232305087831 4.424420050817945 ...
      4.6028652667164938 4.785059376196835 4.9710028287332237 5.1606960308316827 ...
      5.3541393508549282 5.5513331232282894 5.7522776521155432 5.9569732146394792 ...
      6.1654200637104273 6.3776184305161889 6.5935685267188484 6.8132705463971917 ...
      7.0367246677677571 7.2639310547128586 7.4948898581398673 7.729601217192732 ...
      7.9680652603337263 8.2102821063110767 8.456251865026 8.7059746383109058 ...
      8.9594505206290211 9.2166795997043085 9.47766195708957 9.7423976686795033 ...
      10.010886805174749 10.283129432502157 10.559125612196011 10.838875401744188 ...
      11.122378854902983 11.409636021983703 11.700646950113997 11.995411683476323 ...
      12.293930263525857 12.596202729189882 12.90222911705032 13.212009461511173 ...
      13.525543794952087 13.842832147869567 14.163874549006723 14.488671025472797 ...
      14.817221602853252 15.149526305311332 15.485585155681816 15.825398175557639 ...
      16.168965385369994 16.516286804462492 16.867362451159867 17.222192342831622 ...
      17.580776495951152 17.94311492615061 18.309207648271911 18.679054676414133 ...
      ];

B7 = [2.9255283175492273E-8 0.018847166214359216 0.056302231956028639 0.10517026742413378 ...
      0.15808246295731304 0.21050756914460986 0.26092378037044861 0.30949072660246507 ...
      0.35701394498012384 0.40441240839886172 0.45250858893178492 0.5019732152113956 ...
      0.55333118745993881 0.60698608679639476 0.66324647661366865 0.72234848990297729 ...
      0.78447365253459267 0.84976242316525 0.91832430576427659 0.99024535575836348 ...
      1.065593750531485 1.1444239369837992 1.2267797359866877 1.3126966810663998 ...
      1.40220379267138 1.495324934145444 1.5920798557222147 1.6924850042396538 ...
      1.7965541556782 1.9042989127533878 2.0157290990022747 2.1308530729284674 ...
      2.2496779799890159 2.3722099559343945 2.498454291835269 2.6284155687509045 ...
      2.7620977682015062 2.899504363247309 3.0406383939399029 3.1855025301147486 ...
      3.3340991238787527 3.48643025366893 3.6424977613850857 3.8023032838062925 ...
      3.9658482792697702 4.1331340504072029 4.3041617635873335 4.4789324655965546 ...
      4.6574470979949947 4.83970650950952 5.0257114667633074 5.2154626635913663 ...
      5.4089607291502668 5.6062062349965789 5.8071997012806973 6.0119416021797951 ...
      6.2204323706745308 6.4326724027582882 6.6486620611544911 6.86840167860643 ...
      7.0918915607947239 7.3191319889296826 7.5501232220592289 7.7848654991274175 ...
      8.0233590408137943 8.2656040511798032 8.5116007191449885 8.7613492198127769 ...
      9.0148497156630416 9.2721023576265029 9.5331072860541717 9.797864631593324 ...
      10.06637451598016 10.33863705275807 10.6146523479294 10.894420500547561 ...
      11.177941603255741 11.465215742777559 11.756243000364531 12.051023452204639 ...
      12.349557169795773 12.651844220287503 12.957884666794167 13.267678568682021 ...
      13.581225981832823 13.898526958886142 14.219581549462184 14.54438980036702 ...
      14.872951755781703 15.205267457436809 15.541336944773526 15.881160255092574 ...
      16.224737423691952 16.572068483994471 16.923153467665909 17.277992404724575 ...
      17.636585323643011 17.998932251442458 18.365033213780613 18.734888235033306 ...
      ];

B8 = [2.5599371206584062E-8 0.017474422048261291 0.052326626878248433 0.099364942912720022 ...
      0.15297471852633665 0.20852954456706743 0.26359576816973279 0.317422677983887 ...
      0.37019465837132987 0.42249238800526684 0.47499415896895036 0.52833867492292452 ...
      0.58307565752410684 0.63965793386213665 0.69845038840209384 0.75974403417839409 ...
      0.82377014955632422 0.89071262879357982 0.96071811615112679 1.0339040628499827 ...
      1.1103650301430978 1.1901775831422492 1.2734040810116387 1.3600956142449423 ...
      1.4502942868408304 1.544034996258504 1.6413468280439545 1.7422541540612662 ...
      1.8467775019051713 1.9549342468944475 2.0667391658403504 2.1822048825824454 ...
      2.3013422283382541 2.4241605346568811 2.5506678727754282 2.6808712501341261 ...
      2.8147767724761166 2.9523897781654203 3.0937149499713481 3.2387564084917142 ...
      3.3875177905476965 3.5400023152249189 3.696212839716829 3.8561519067160241 ...
      4.0198217847730415 4.1872245027815405 4.3583618795400207 4.5332355491719367 ...
      4.7118469830500764 4.8941975087607128 5.080288326552977 5.2701205236454109 ...
      5.4636950867011667 5.6610129127336339 5.8620748186630758 6.0668815497107795 ...
      6.2754337867888808 6.4877321530202448 6.7037772195030607 6.9235695104181154 ...
      7.1471095075626847 7.3743976543832073 7.6054343595688287 7.8402200002594649 ...
      8.07875492491475 8.3210394558841134 8.567073891712921 8.81685850921514 ...
      9.0703935653390442 9.3276792988492012 9.58871593184509 9.8535036711341384 ...
      10.122042709474913 10.394333226704266 10.670375390760617 10.950169358614154 ...
      11.233715277113452 11.521013283757 11.812063507397138 12.106866068883033 ...
      12.405421081648669 12.707728652251172 13.013788880864141 13.32360186173028 ...
      13.63716768357706 13.954486429998903 14.27555817980881 14.600383007362314 ...
      14.928960982856072 15.261292172603531 15.597376639289443 15.937214442205198 ...
      16.280805637466575 16.628150278215383 16.979248414806346 17.334100094980403 ...
      17.692705364025674 18.055064264926909 18.421176838504483 18.791043123543677 ...
      ];

B9 = [2.280064796632299E-8 0.016362982921077733 0.049097816135343912 0.094017747658539055 ...
      0.14710212926762414 0.20415606613288764 0.26236903947605772 0.32034844609230145 ...
      0.37769940123803458 0.43459467480237657 0.4914690966156321 0.54883882769220316 ...
      0.60720846972418085 0.66703059775342177 0.72869333029312988 0.79252142967603445 ...
      0.85878304680022688 0.9276981291544476 0.9994466459419844 1.0741758869236178 ...
      1.1520066296606448 1.2330382125617092 1.3173526456506377 1.4050179147349455 ...
      1.4960906271965477 1.5906181284418168 1.688640196298244 1.7901904002904843 ...
      1.8952971951955035 2.0039848038154133 2.1162739322700244 2.2321823518820736 ...
      2.3517253744690074 2.4749162421711755 2.6017664485036986 2.7322860038475967 ...
      2.8664836558797249 3.0043670733121379 3.1459429996385477 3.2912173822669675 ...
      3.4401954813750528 3.5928819619974282 3.7492809721956406 3.90939620963498 ...
      4.0732309784702956 4.2407882381030326 4.4120706450972191 4.587080589319446 ...
      4.7658202251868049 4.9482914987588007 5.1344961712880943 5.3244358397452984 ...
      5.5181119547509105 5.7155258362794727 5.91667868744463 6.1215716066268024 ...
      6.3302055981659935 6.542581581809328 6.7587004010754637 6.9785628306747309 ...
      7.2021695831042374 7.4295213145206285 7.6606186299790808 7.895462088115182 ...
      8.1340522053360189 8.3763894595782133 8.6224742936830339 8.8723071184324063 ...
      9.1258883152840138 9.38321823883897 9.6442972190714684 9.9091255633461532 ...
      10.177703558245913 10.450031471230163 10.726109552141274 11.005938034574742 ...
      11.289517137127032 11.576847064533299 11.867928008706034 12.162760149684221 ...
      12.461343656501773 12.763678687982955 13.069765393471682 13.379603913500933 ...
      13.693194380407759 14.010536918898946 14.331631646571715 14.656478674393551 ...
      14.985078107144719 15.317430043826853 15.65353457804037 15.993391798333567 ...
      16.337001788525711 16.684364628006307 17.03548039201258 17.390349151886877 ...
      17.748970975315721 18.111345926551934 18.477474066621202 18.847355453514293 ...
      ];

D = 4.6906779762565158E-8;

E1 = 2.5E+9;

E2 = 2.78E+10;

E3 = 2.45E+9;

F = [10 10110.90909090909 20211.81818181818 30312.727272727272 40413.63636363636 ...
     50514.545454545456 60615.454545454544 70716.363636363632 80817.272727272721 ...
     90918.181818181823 101019.09090909091 111120 121220.90909090909 131321.81818181818 ...
     141422.72727272726 151523.63636363635 161624.54545454544 171725.45454545456 ...
     181826.36363636365 191927.27272727274 202028.18181818182 212129.09090909091 ...
     222230 232330.90909090909 242431.81818181818 252532.72727272726 262633.63636363635 ...
     272734.54545454547 282835.45454545453 292936.36363636365 303037.27272727271 ...
     313138.18181818182 323239.09090909088 333340 343440.90909090912 353541.81818181818 ...
     363642.72727272729 373743.63636363635 383844.54545454547 393945.45454545453 ...
     404046.36363636365 414147.27272727271 424248.18181818182 434349.09090909088 ...
     444450 454550.90909090912 464651.81818181818 474752.72727272729 484853.63636363635 ...
     494954.54545454547 505055.45454545453 515156.36363636365 525257.27272727271 ...
     535358.18181818177 545459.09090909094 555560 565660.90909090906 575761.81818181823 ...
     585862.72727272729 595963.63636363635 606064.54545454541 616165.45454545459 ...
     626266.36363636365 636367.27272727271 646468.18181818177 656569.09090909094 ...
     666670 676770.90909090906 686871.81818181823 696972.72727272729 707073.63636363635 ...
     717174.54545454541 727275.45454545459 737376.36363636365 747477.27272727271 ...
     757578.18181818177 767679.09090909094 777780 787880.90909090906 797981.81818181823 ...
     808082.72727272729 818183.63636363635 828284.54545454541 838385.45454545459 ...
     848486.36363636365 858587.27272727271 868688.18181818177 878789.09090909094 ...
     888890 898990.90909090906 909091.81818181823 919192.72727272729 929293.63636363635 ...
     939394.54545454541 949495.45454545459 959596.36363636365 969697.27272727271 ...
     979798.18181818177 989899.09090909094 1.0E+6];

Fr00 = 24;

Fr010 = 3780.06067615805;

Fr0100 = 80373.924748103862;

Fr020 = 10545.429753925788;

Fr030 = 18397.440084666694;

Fr040 = 26763.725524520923;

Fr050 = 35413.859616810849;

Fr060 = 44237.187620112294;

Fr070 = 53173.956738159912;

Fr080 = 62189.0658844067;

Fr090 = 71260.541368135513;

FrN0 = 9;

FrN10 = 1049.80677635431;

FrN100 = 10417.0677635431;

FrN20 = 2090.61355270862;

FrN30 = 3131.4203290629303;

FrN40 = 4172.22710541724;

FrN50 = 5213.03388177155;

FrN60 = 6253.8406581258605;

FrN70 = 7294.6474344801709;

FrN80 = 8335.45421083448;

FrN90 = 9376.2609871887926;

Ps = 1;

Ps0 = 1;

Psat = 0.023060749597593345;

T = 293.15;

T0 = 293.15;

T01 = 273.16;

Z0 = 419;

a = 0.00013;

alpha0 = saveVarsMat.alpha0; % <1x1 function_handle> unsupported class

alpha1 = saveVarsMat.alpha1; % <1x1 function_handle> unsupported class

alpha10 = saveVarsMat.alpha10; % <1x1 function_handle> unsupported class

alpha2 = saveVarsMat.alpha2; % <1x1 function_handle> unsupported class

alpha3 = saveVarsMat.alpha3; % <1x1 function_handle> unsupported class

alpha4 = saveVarsMat.alpha4; % <1x1 function_handle> unsupported class

alpha5 = saveVarsMat.alpha5; % <1x1 function_handle> unsupported class

alpha6 = saveVarsMat.alpha6; % <1x1 function_handle> unsupported class

alpha7 = saveVarsMat.alpha7; % <1x1 function_handle> unsupported class

alpha8 = saveVarsMat.alpha8; % <1x1 function_handle> unsupported class

alpha9 = saveVarsMat.alpha9; % <1x1 function_handle> unsupported class

c = 343;

co = 343;

dist1 = -6.2135978267359222E-6;

dist2 = -6.0135978267359224E-6;

dist3 = -5.5135978267359221E-6;

h0 = 0;

h1 = 0.23060749597593344;

h10 = 2.3060749597593344;

h2 = 0.46121499195186688;

h3 = 0.69182248792780032;

h4 = 0.92242998390373376;

h5 = 1.1530374798796672;

h6 = 1.3836449758556006;

h7 = 1.6142524718315341;

h8 = 1.8448599678074675;

h9 = 2.0754674637834012;

hh1 = 1.5E-6;

hh2 = 1.6999999999999998E-6;

hh3 = 2.2E-6;

hr = [0 10 20 30 40 50 60 70 80 90 100];

lambda = 0.0005346193987233117;

lambda00 = 3.1968734726291563;

mueff = 0.041673999999999996;

plate1 = 2.8267752148349166E+9;

plate2 = 3.0352658587182007E+10;

plate3 = 2.770239710538218E+9;

r = 0.1;

ro1 = 1420;

ro2 = 2700;

ro3 = 1780;

t1 = 2E-7;

t2 = 1.4999999999999999E-5;

t3 = 5E-7;

v1 = 0.34;

v2 = 0.29;

v3 = 0.34;

wn = 641577.91658719268;

z1 = 1E-7;

z2 = 7.6999999999999991E-6;

z3 = 1.545E-5;

zdist = 7.7364021732640772E-6;

zp = 7.7135978267359222E-6;

clear saveVarsMat;