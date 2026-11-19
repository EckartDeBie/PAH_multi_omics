#Multi-omics in PAH analyses scripts 

#Note: some parts of the code have been anonymised (patient IDs removed, exact file path names removed)

############################################################################
#file preparation
#MOFA with missing CIBER adjusted and AABS
#EdB started here 01/11/24
library(tidyverse)
library(survival)
library(survminer)
library(MOFA2)
library(fgsea)
library(msigdbr)
library(data.table)
library(pheatmap)
library(Publish)
getwd()

#import the data
#metbolomics
metabolomics1 <- read_rds("mounted-data-readonly/PAH-omics/omics/Metabolomics/metablomics_normalisd_z_score_metabolomics.rds")
#big problem! Many vars were renamed --> correct this!

metabolomics1 <- metabolomics1 %>%
  rename( `X...11429` 	=	 `X...58`	,
          `X...11564` 	=	 `X...59`	,
          `X...12127` 	=	 `X...60`	,
          `X...12472` 	=	 `X...61`	,
          `X...12688` 	=	 `X...62`	,
          `X...12739` 	=	 `X...63`	,
          `X...12824` 	=	 `X...64`	,
          `X...13737` 	=	 `X...65`	,
          `X...15503` 	=	 `X...66`	,
          `X...17327` 	=	 `X...67`	,
          `X...21796` 	=	 `X...68`	,
          `X...21829` 	=	 `X...69`	,
          `X...23765` 	=	 `X...70`	,
          `X...24020` 	=	 `X...71`	,
          `X...24041` 	=	 `X...72`	,
          `X...24307` 	=	 `X...73`	,
          `X...24411` 	=	 `X...74`	,
          `X...24513` 	=	 `X...75`	,
          `X...24527` 	=	 `X...76`	,
          `X...24678` 	=	 `X...77`	,
          `X...24728` 	=	 `X...78`	,
          `X...24766` 	=	 `X...79`	,
          `X...01911` 	=	 `X...973`	,
          `X...02249` 	=	 `X...974`	,
          `X...02269` 	=	 `X...975`	,
          `X...07765` 	=	 `X...976`	,
          `X...09789` 	=	 `X...977`	,
          `X...10458` 	=	 `X...978`	,
          `X...11261` 	=	 `X...979`	,
          `X...11299` 	=	 `X...980`	,
          `X...11305` 	=	 `X...981`	,
          `X...11308` 	=	 `X...982`	,
          `X...11315` 	=	 `X...983`	,
          `X...11334` 	=	 `X...984`	,
          `X...11357` 	=	 `X...985`	,
          `X...11372` 	=	 `X...986`	,
          `X...11378` 	=	 `X...987`	,
          `X...11381` 	=	 `X...988`	,
          `X...11407` 	=	 `X...989`	,
          `X...11438` 	=	 `X...990`	,
          `X...11440` 	=	 `X...991`	,
          `X...11441` 	=	 `X...992`	,
          `X...11442` 	=	 `X...993`	,
          `X...11444` 	=	 `X...994`	,
          `X...11452` 	=	 `X...995`	,
          `X...11470` 	=	 `X...996`	,
          `X...11478` 	=	 `X...997`	,
          `X...11483` 	=	 `X...998`	,
          `X...11485` 	=	 `X...999`	,
          `X...11491` 	=	 `X...1000`	,
          `X...11522` 	=	 `X...1001`	,
          `X...11530` 	=	 `X...1002`	,
          `X...11537` 	=	 `X...1003`	,
          `X...11538` 	=	 `X...1004`	,
          `X...11540` 	=	 `X...1005`	,
          `X...11549` 	=	 `X...1006`	,
          `X...11632` 	=	 `X...1007`	,
          `X...11787` 	=	 `X...1008`	,
          `X...11795` 	=	 `X...1009`	,
          `X...11838` 	=	 `X...1010`	,
          `X...11843` 	=	 `X...1011`	,
          `X...11847` 	=	 `X...1012`	,
          `X...11849` 	=	 `X...1013`	,
          `X...11850` 	=	 `X...1014`	,
          `X...11852` 	=	 `X...1015`	,
          `X...11858` 	=	 `X...1016`	,
          `X...11861` 	=	 `X...1017`	,
          `X...11871` 	=	 `X...1018`	,
          `X...11880` 	=	 `X...1019`	,
          `X...11905` 	=	 `X...1020`	,
          `X...12007` 	=	 `X...1021`	,
          `X...12013` 	=	 `X...1022`	,
          `X...12015` 	=	 `X...1023`	,
          `X...12026` 	=	 `X...1024`	,
          `X...12063` 	=	 `X...1025`	,
          `X...12100` 	=	 `X...1026`	,
          `X...12101` 	=	 `X...1027`	,
          `X...12104` 	=	 `X...1028`	,
          `X...12111` 	=	 `X...1029`	,
          `X...12117` 	=	 `X...1030`	,
          `X...12123` 	=	 `X...1031`	,
          `X...12125` 	=	 `X...1032`	,
          `X...12126` 	=	 `X...1033`	,
          `X...12170` 	=	 `X...1034`	,
          `X...12193` 	=	 `X...1035`	,
          `X...12199` 	=	 `X...1036`	,
          `X...12206` 	=	 `X...1037`	,
          `X...12212` 	=	 `X...1038`	,
          `X...12216` 	=	 `X...1039`	,
          `X...12221` 	=	 `X...1040`	,
          `X...12230` 	=	 `X...1041`	,
          `X...12231` 	=	 `X...1042`	,
          `X...12257` 	=	 `X...1043`	,
          `X...12261` 	=	 `X...1044`	,
          `X...12263` 	=	 `X...1045`	,
          `X...12283` 	=	 `X...1046`,
          `X...12306` 	=	 `X...1047`,
          `X...12329` 	=	 `X...1048`	,
          `X...12339` 	=	 `X...1049`	,
          `X...12407` 	=	 `X...1050`	,
          `X...12410` 	=	 `X...1051`	,
          `X...12411` 	=	 `X...1052`	,
          `X...12442` 	=	 `X...1053`	,
          `X...12450` 	=	 `X...1054`	,
          `X...12456` 	=	 `X...1055`	,
          `X...12459` 	=	 `X...1056`	,
          `X...12462` 	=	 `X...1057`	,
          `X...12524` 	=	 `X...1058`	,
          `X...12543` 	=	 `X...1059`	,
          `X...12544` 	=	 `X...1060`	,
          `X...12565` 	=	 `X...1061`	,
          `X...12636` 	=	 `X...1062`	,
          `X...12680` 	=	 `X...1063`	,
          `X...12681` 	=	 `X...1064`	,
          `X...12686` 	=	 `X...1065`	,
          `X...12701` 	=	 `X...1066`	,
          `X...12707` 	=	 `X...1067`	,
          `X...12708` 	=	 `X...1068`	,
          `X...12712` 	=	 `X...1069`	,
          `X...12714` 	=	 `X...1070`	,
          `X...12715` 	=	 `X...1071`	,
          `X...12718` 	=	 `X...1072`	,
          `X...12726` 	=	 `X...1073`	,
          `X...12729` 	=	 `X...1074`	,
          `X...12730` 	=	 `X...1075`	,
          `X...12731` 	=	 `X...1076`	,
          `X...12733` 	=	 `X...1077`	,
          `X...12738` 	=	 `X...1078`	,
          `X...12740` 	=	 `X...1079`	,
          `X...12748` 	=	 `X...1080`	,
          `X...12753` 	=	 `X...1081`	,
          `X...12798` 	=	 `X...1082`	,
          `X...12812` 	=	 `X...1083`	,
          `X...12814` 	=	 `X...1084`	,
          `X...12816` 	=	 `X...1085`	,
          `X...12818` 	=	 `X...1086`	,
          `X...12819` 	=	 `X...1087`	,
          `X...12822` 	=	 `X...1088`	,
          `X...12830` 	=	 `X...1089`	,
          `X...12844` 	=	 `X...1090`	,
          `X...12846` 	=	 `X...1091`	,
          `X...12847` 	=	 `X...1092`	,
          `X...12849` 	=	 `X...1093`	,
          `X...12851` 	=	 `X...1094`	,
          `X...12879` 	=	 `X...1095`	,
          `X...12906` 	=	 `X...1096`	,
          `X...13007` 	=	 `X...1097`	,
          `X...13431` 	=	 `X...1098`	,
          `X...13435` 	=	 `X...1099`	,
          `X...13529` 	=	 `X...1100`	,
          `X...13553` 	=	 `X...1101`	,
          `X...13658` 	=	 `X...1102`	,
          `X...13684` 	=	 `X...1103`	,
          `X...13688` 	=	 `X...1104`	,
          `X...13689` 	=	 `X...1105`	,
          `X...13703` 	=	 `X...1106`	,
          `X...13722` 	=	 `X...1107`	,
          `X...13723` 	=	 `X...1108`	,
          `X...13726` 	=	 `X...1109`	,
          `X...13729` 	=	 `X...1110`	,
          `X...13835` 	=	 `X...1111`	,
          `X...13837` 	=	 `X...1112`	,
          `X...13844` 	=	 `X...1113`	,
          `X...13846` 	=	 `X...1114`	,
          `X...13866` 	=	 `X...1115`	,
          `X...13891` 	=	 `X...1116`	,
          `X...14056` 	=	 `X...1117`	,
          `X...14095` 	=	 `X...1118`	,
          `X...14113` 	=	 `X...1119`	,
          `X...14141` 	=	 `X...1120`	,
          `X...14302` 	=	 `X...1121`	,
          `X...14314` 	=	 `X...1122`	,
          `X...14364` 	=	 `X...1123`	,
          `X...14568` 	=	 `X...1124`	,
          `X...14626` 	=	 `X...1125`	,
          `X...14658` 	=	 `X...1126`	,
          `X...14662` 	=	 `X...1127`	,
          `X...14838` 	=	 `X...1128`	,
          `X...14904` 	=	 `X...1129`	,
          `X...14939` 	=	 `X...1130`	,
          `X...15220` 	=	 `X...1131`	,
          `X...15245` 	=	 `X...1132`	,
          `X...15469` 	=	 `X...1133`	,
          `X...15486` 	=	 `X...1134`	,
          `X...15492` 	=	 `X...1135`	,
          `X...15497` 	=	 `X...1136`	,
          `X...15666` 	=	 `X...1137`	,
          `X...16071` 	=	 `X...1138`	,
          `X...16083` 	=	 `X...1139`	,
          `X...16124` 	=	 `X...1140`	,
          `X...16397` 	=	 `X...1141`	,
          `X...16480` 	=	 `X...1142`	,
          `X...16570` 	=	 `X...1143`	,
          `X...16576` 	=	 `X...1144`	,
          `X...16580` 	=	 `X...1145`	,
          `X...16649` 	=	 `X...1146`	,
          `X...16654` 	=	 `X...1147`	,
          `X...16935` 	=	 `X...1148`	,
          `X...16938` 	=	 `X...1149`	,
          `X...16944` 	=	 `X...1150`	,
          `X...16946` 	=	 `X...1151`	,
          `X...16947` 	=	 `X...1152`	,
          `X...16964` 	=	 `X...1153`	,
          `X...17010` 	=	 `X...1154`	,
          `X...17145` 	=	 `X...1155`	,
          `X...17146` 	=	 `X...1156`	,
          `X...17166` 	=	 `X...1157`	,
          `X...17185` 	=	 `X...1158`	,
          `X...17189` 	=	 `X...1159`	,
          `X...17269` 	=	 `X...1160`	,
          `X...17299` 	=	 `X...1161`	,
          `X...17325` 	=	 `X...1162`	,
          `X...17328` 	=	 `X...1163`	,
          `X...17335` 	=	 `X...1164`	,
          `X...17337` 	=	 `X...1165`	,
          `X...17340` 	=	 `X...1166`	,
          `X...17343` 	=	 `X...1167`	,
          `X...17346` 	=	 `X...1168`	,
          `X...17348` 	=	 `X...1169`	,
          `X...17350` 	=	 `X...1170`	,
          `X...17351` 	=	 `X...1171`	,
          `X...17353` 	=	 `X...1172`	,
          `X...17354` 	=	 `X...1173`	,
          `X...17355` 	=	 `X...1174`	,
          `X...17359` 	=	 `X...1175`	,
          `X...17361` 	=	 `X...1176`	,
          `X...17367` 	=	 `X...1177`	,
          `X...17393` 	=	 `X...1178`	,
          `X...17398` 	=	 `X...1179`	,
          `X...17438` 	=	 `X...1180`	,
          `X...17469` 	=	 `X...1181`	,
          `X...17653` 	=	 `X...1182`	,
          `X...17654` 	=	 `X...1183`	,
          `X...17655` 	=	 `X...1184`	,
          `X...17673` 	=	 `X...1185`	,
          `X...17676` 	=	 `X...1186`	,
          `X...17677` 	=	 `X...1187`	,
          `X...17685` 	=	 `X...1188`	,
          `X...17686` 	=	 `X...1189`	,
          `X...17690` 	=	 `X...1190`	,
          `X...17692` 	=	 `X...1191`	,
          `X...17694` 	=	 `X...1192`	,
          `X...17704` 	=	 `X...1193`	,
          `X...17705` 	=	 `X...1194`	,
          `X...17765` 	=	 `X...1195`	,
          `X...17969` 	=	 `X...1196`	,
          `X...18165` 	=	 `X...1197`	,
          `X...18249` 	=	 `X...1198`	,
          `X...18345` 	=	 `X...1199`	,
          `X...18603` 	=	 `X...1200`	,
          `X...18779` 	=	 `X...1201`	,
          `X...18886` 	=	 `X...1202`	,
          `X...18887` 	=	 `X...1203`	,
          `X...18889` 	=	 `X...1204`	,
          `X...18899` 	=	 `X...1205`	,
          `X...18901` 	=	 `X...1206`	,
          `X...18913` 	=	 `X...1207`	,
          `X...18914` 	=	 `X...1208`	,
          `X...18921` 	=	 `X...1209`	,
          `X...18922` 	=	 `X...1210`	,
          `X...19141` 	=	 `X...1211`	,
          `X...19183` 	=	 `X...1212`	,
          `X...19299` 	=	 `X...1213`	,
          `X...19434` 	=	 `X...1214`	,
          `X...19438` 	=	 `X...1215`	,
          `X...19455` 	=	 `X...1216`	,
          `X...19561` 	=	 `X...1217`	,
          `X...20529` 	=	 `X...1218`	,
          `X...21258` 	=	 `X...1219`	,
          `X...21285` 	=	 `X...1220`	,
          `X...21286` 	=	 `X...1221`	,
          `X...21295` 	=	 `X...1222`	,
          `X...21310` 	=	 `X...1223`	,
          `X...21312` 	=	 `X...1224`	,
          `X...21315` 	=	 `X...1225`	,
          `X...21318` 	=	 `X...1226`	,
          `X...21319` 	=	 `X...1227`	,
          `X...21339` 	=	 `X...1228`	,
          `X...21341` 	=	 `X...1229`	,
          `X...21343` 	=	 `X...1230`	,
          `X...21353` 	=	 `X...1231`	,
          `X...21358` 	=	 `X...1232`	,
          `X...21364` 	=	 `X...1233`	,
          `X...21365` 	=	 `X...1234`	,
          `X...21383` 	=	 `X...1235`	,
          `X...21410` 	=	 `X...1236`	,
          `X...21411` 	=	 `X...1237`	,
          `X...21441` 	=	 `X...1238`	,
          `X...21442` 	=	 `X...1239`	,
          `X...21444` 	=	 `X...1240`	,
          `X...21448` 	=	 `X...1241`	,
          `X...21470` 	=	 `X...1242`	,
          `X...21607` 	=	 `X...1243`	,
          `X...21626` 	=	 `X...1244`	,
          `X...21657` 	=	 `X...1245`	,
          `X...21658` 	=	 `X...1246`	,
          `X...21661` 	=	 `X...1247`	,
          `X...21662` 	=	 `X...1248`	,
          `X...21666` 	=	 `X...1249`	,
          `X...21668` 	=	 `X...1250`	,
          `X...21729` 	=	 `X...1251`	,
          `X...21733` 	=	 `X...1252`	,
          `X...21735` 	=	 `X...1253`	,
          `X...21736` 	=	 `X...1254`	,
          `X...21740` 	=	 `X...1255`	,
          `X...21752` 	=	 `X...1256`	,
          `X...21792` 	=	 `X...1257`	,
          `X...21795` 	=	 `X...1258`	,
          `X...21803` 	=	 `X...1259`	,
          `X...21807` 	=	 `X...1260`	,
          `X...21815` 	=	 `X...1261`	,
          `X...21816` 	=	 `X...1262`	,
          `X...21821` 	=	 `X...1263`	,
          `X...21831` 	=	 `X...1264`	,
          `X...21834` 	=	 `X...1265`	,
          `X...21838` 	=	 `X...1266`	,
          `X...21849` 	=	 `X...1267`	,
          `X...21935` 	=	 `X...1268`	,
          `X...21959` 	=	 `X...1269`	,
          `X...22035` 	=	 `X...1270`	,
          `X...22102` 	=	 `X...1271`	,
          `X...22143` 	=	 `X...1272`	,
          `X...22145` 	=	 `X...1273`	,
          `X...22147` 	=	 `X...1274`	,
          `X...22162` 	=	 `X...1275`	,
          `X...22379` 	=	 `X...1276`	,
          `X...22475` 	=	 `X...1277`	,
          `X...22508` 	=	 `X...1278`	,
          `X...22509` 	=	 `X...1279`	,
          `X...22515` 	=	 `X...1280`	,
          `X...22520` 	=	 `X...1281`	,
          `X...22716` 	=	 `X...1282`	,
          `X...22764` 	=	 `X...1283`	,
          `X...22771` 	=	 `X...1284`	,
          `X...22775` 	=	 `X...1285`	,
          `X...22776` 	=	 `X...1286`	,
          `X...22789` 	=	 `X...1287`	,
          `X...22834` 	=	 `X...1288`	,
          `X...22836` 	=	 `X...1289`	,
          `X...22850` 	=	 `X...1290`	,
          `X...23145` 	=	 `X...1291`	,
          `X...23159` 	=	 `X...1292`	,
          `X...23196` 	=	 `X...1293`	,
          `X...23276` 	=	 `X...1294`	,
          `X...23283` 	=	 `X...1295`	,
          `X...23291` 	=	 `X...1296`	,
          `X...23293` 	=	 `X...1297`	,
          `X...23294` 	=	 `X...1298`	,
          `X...23295` 	=	 `X...1299`	,
          `X...23297` 	=	 `X...1300`	,
          `X...23314` 	=	 `X...1301`	,
          `X...23369` 	=	 `X...1302`	,
          `X...23583` 	=	 `X...1303`	,
          `X...23585` 	=	 `X...1304`	,
          `X...23587` 	=	 `X...1305`	,
          `X...23588` 	=	 `X...1306`	,
          `X...23590` 	=	 `X...1307`	,
          `X...23636` 	=	 `X...1308`	,
          `X...23637` 	=	 `X...1309`	,
          `X...23639` 	=	 `X...1310`	,
          `X...23641` 	=	 `X...1311`	,
          `X...23644` 	=	 `X...1312`	,
          `X...23648` 	=	 `X...1313`	,
          `X...23649` 	=	 `X...1314`	,
          `X...23652` 	=	 `X...1315`	,
          `X...23655` 	=	 `X...1316`	,
          `X...23659` 	=	 `X...1317`	,
          `X...23669` 	=	 `X...1318`	,
          `X...23680` 	=	 `X...1319`	,
          `X...23705` 	=	 `X...1320`	,
          `X...23729` 	=	 `X...1321`	,
          `X...23739` 	=	 `X...1322`	,
          `X...23749` 	=	 `X...1323`	,
          `X...23756` 	=	 `X...1324`	,
          `X...23776` 	=	 `X...1325`	,
          `X...23780` 	=	 `X...1326`	,
          `X...23782` 	=	 `X...1327`	,
          `X...23787` 	=	 `X...1328`	,
          `X...23868` 	=	 `X...1329`	,
          `X...23984` 	=	 `X...1330`	,
          `X...23997` 	=	 `X...1331`	,
          `X...24023` 	=	 `X...1332`	,
          `X...24024` 	=	 `X...1333`	,
          `X...24027` 	=	 `X...1334`	,
          `X...24049` 	=	 `X...1335`	,
          `X...24054` 	=	 `X...1336`	,
          `X...24061` 	=	 `X...1337`	,
          `X...24097` 	=	 `X...1338`	,
          `X...24125` 	=	 `X...1339`	,
          `X...24204` 	=	 `X...1340`	,
          `X...24220` 	=	 `X...1341`	,
          `X...24240` 	=	 `X...1342`	,
          `X...24243` 	=	 `X...1343`	,
          `X...24260` 	=	 `X...1344`	,
          `X...24261` 	=	 `X...1345`	,
          `X...24263` 	=	 `X...1346`	,
          `X...24278` 	=	 `X...1347`	,
          `X...24293` 	=	 `X...1348`	,
          `X...24295` 	=	 `X...1349`	,
          `X...24306` 	=	 `X...1350`	,
          `X...24309` 	=	 `X...1351`	,
          `X...24329` 	=	 `X...1352`	,
          `X...24334` 	=	 `X...1353`	,
          `X...24337` 	=	 `X...1354`	,
          `X...24348` 	=	 `X...1355`	,
          `X...24352` 	=	 `X...1356`	,
          `X...24418` 	=	 `X...1357`	,
          `X...24422` 	=	 `X...1358`	,
          `X...24425` 	=	 `X...1359`	,
          `X...24431` 	=	 `X...1360`	,
          `X...24432` 	=	 `X...1361`	,
          `X...24438` 	=	 `X...1362`	,
          `X...24439` 	=	 `X...1363`	,
          `X...24449` 	=	 `X...1364`	,
          `X...24452` 	=	 `X...1365`	,
          `X...24455` 	=	 `X...1366`	,
          `X...24456` 	=	 `X...1367`	,
          `X...24475` 	=	 `X...1368`	,
          `X...24512` 	=	 `X...1369`	,
          `X...24514` 	=	 `X...1370`	,
          `X...24518` 	=	 `X...1371`	,
          `X...24539` 	=	 `X...1372`	,
          `X...24542` 	=	 `X...1373`	,
          `X...24544` 	=	 `X...1374`	,
          `X...24546` 	=	 `X...1375`	,
          `X...24549` 	=	 `X...1376`	,
          `X...24550` 	=	 `X...1377`	,
          `X...24551` 	=	 `X...1378`	,
          `X...24553` 	=	 `X...1379`	,
          `X...24554` 	=	 `X...1380`	,
          `X...24557` 	=	 `X...1381`	,
          `X...24558` 	=	 `X...1382`	,
          `X...24559` 	=	 `X...1383`	,
          `X...24560` 	=	 `X...1384`	,
          `X...24565` 	=	 `X...1385`	,
          `X...24568` 	=	 `X...1386`	,
          `X...24571` 	=	 `X...1387`	,
          `X...24580` 	=	 `X...1388`	,
          `X...24585` 	=	 `X...1389`	,
          `X...24588` 	=	 `X...1390`	,
          `X...24625` 	=	 `X...1391`	,
          `X...24637` 	=	 `X...1392`	,
          `X...24655` 	=	 `X...1393`	,
          `X...24657` 	=	 `X...1394`	,
          `X...24669` 	=	 `X...1395`	,
          `X...24699` 	=	 `X...1396`	,
          `X...24736` 	=	 `X...1397`	,
          `X...24738` 	=	 `X...1398`	,
          `X...24747` 	=	 `X...1399`	,
          `X...24748` 	=	 `X...1400`	,
          `X...24749` 	=	 `X...1401`	,
          `X...24757` 	=	 `X...1402`	,
          `X...24760` 	=	 `X...1403`	,
          `X...24761` 	=	 `X...1404`	,
          `X...24762` 	=	 `X...1405`	,
          `X...24765` 	=	 `X...1406`	,
          `X...24803` 	=	 `X...1407`	,
          `X...24804` 	=	 `X...1408`	,
          `X...24806` 	=	 `X...1409`	,
          `X...24809` 	=	 `X...1410`	,
          `X...24811` 	=	 `X...1411`	,
          `X...24812` 	=	 `X...1412`	,
          `X...24813` 	=	 `X...1413`	,
          `X...24831` 	=	 `X...1414`	,
          `X...24832` 	=	 `X...1415`	,
          `X...24845` 	=	 `X...1416`	,
          `X...24849` 	=	 `X...1417`	)


#RNA
residuals_ciber_corrected <- readRDS("/mnt/session_data/residuals_ciber_corrected.rds")
rna <- residuals_ciber_corrected

#Z-scored residuals for proteomics)
proteomics_correct <- readRDS("~/session_data/mounted-data-readonly/PAH-omics/omics/Proteomics/proteomics_correct.rds")

#AABs
Autoantibody_levels_patients_log_transformed <- readRDS("/mnt/session_data/mounted-data-readonly/PAH-omics/omics/Autoantibodies/Autoantibody_levels_patients_log_transformed.rds")

#WBCs
PCs_WBCs_RNA <- readRDS("/mnt/session_data/mounted-data-readonly/PCs_WBCs_RNA.rds")
#=============================================================================================================
#get IDs sorted
pheno_rna_seq <- readRDS("/mnt/session_data/mounted-data-readonly/PAH-omics/omics/RNA_seq/pheno_rna_seq.rds")
id <- pheno_rna_seq %>% select(SampleID, OpenClinicaID)
follow_up_data <- readRDS("~/session_data/mounted-data-readonly/PAH-OpenClinica/20240726-freeze/follow_up_data.rds")
id2 <- follow_up_data %>% select(id, id_cohort)
id <- merge(id2, id, by.x='id', by.y='OpenClinicaID', all=T)
id <- id %>% filter(!is.na(id_cohort))
id <- unique(id)
#============================================================================================================
#format proteomics
proteomics_correct[proteomics_correct == ""] <- NA
colnames(proteomics_correct)[is.na(colnames(proteomics_correct)) | colnames(proteomics_correct) == ""] <- "Unnamed_Column"
proteomics_correct <- proteomics_correct %>% filter(!is.na(OpenClinica.ID))

#select relevant proteomics
prot1 <- proteomics_correct %>% filter(Event == 'PATIENT_1')
prot1 <- prot1[,-c(1,3:64)]
prot1 <- t(prot1)
colnames(prot1) <- prot1[1,]
prot1 <- prot1[-1,]
#NB: including FC.mouse here.....
#need to check if appropriate!
prot1 <- as.data.frame(prot1)
prot1[,1:463] <- lapply(prot1[,1:463], as.numeric)
prot1 <- as.data.frame(prot1)
#this is complete

#now metabolomics
#NOTE: there is a substantial number of unnamed metabolites
metabolomics <- merge(id2, metabolomics1, by.x='id_cohort', by.y='SUBJECT_ID')
metab1 <- metabolomics %>% select(id, where(~ !is.character(.)))
metab1 <- unique(metab1)
#the number used to correspond to exact number in Rhodes' paper, so should be correct!
metab1 <- t(metab1)
colnames(metab1) <- metab1[1,]
metab1 <- metab1[-1,]
metab1 <- as.data.frame(metab1)
metab1[,1:288] <- lapply(metab1[,1:288], as.numeric)
#quite a lot of vars are NA --> remove for now
#metab1 <- na.omit(metab1)
#metab1 <- as.data.frame(metab1)
#metabolomics are now in shape

#check NAs
#this has worked!! Now we just need to filter on score
row_na_sums <- apply(metab1, 1, function(x) sum(is.na(x)))
metab1$na_count <- row_na_sums

#NA count here is 0?
metab1 <-metab1 %>% select(-na_count)


#now select this per dataset
rna1 <- rna
pheno_rna_seq <- readRDS("/mnt/session_data/mounted-data-readonly/PAH-omics/omics/RNA_seq/pheno_rna_seq.rds")
keep <- pheno_rna_seq %>% filter(PAH==1)

rna1 <- merge(rna1, keep[,c(9,10)], by.x='ID', by.y='CohortID')

#prevent overlapping IDs
rna1 <- rna1 %>% distinct(OpenClinicaID, .keep_all = TRUE)
rna1 <- rna1[,-1]
rna1$OpenClinicaID <- as.character(rna1$OpenClinicaID)
rna1 <- rna1[grepl("^OC", rna1$OpenClinicaID), ]
#has finally worked!

#now transpose and prepare
rna1 <- as.data.frame(t(rna1))
#363 patients retained. 
rna2 <- rna1
colnames(rna2) <- rna2[15083,]
rna2 <- rna2[-c(15083),]
#formatted correctly now
rna2 <- as.data.frame(rna2)
rna2[,1:363] <- lapply(rna2[,1:363], as.numeric)
rna2 <- as.data.frame(rna2)

#include CIBER too
ciber <- PCs_WBCs_RNA[,c(1,17:38)]
ciber <- merge(ciber, id, by.x = 'CohortID', by.y='id_cohort')
ciber <- ciber[,-1]
ciber <- ciber %>% distinct(id, .keep_all = TRUE)
cibter <- ciber %>% filter(id %in% colnames(rna2))
ciber <- cibter[, colSums(cibter != 0) > 0]

c2 <- t(ciber)
colnames(c2) <- c2[19,]
c2 <- c2[-c(19,20),]

c2 <- as.data.frame(c2)
c2[,1:363] <- lapply(c2[,1:363], as.numeric)

#now get aabs in shape
aab1 <- Autoantibody_levels_patients_log_transformed
aab1 <- merge(aab1, id, by.x='sample', by.y='id_cohort')
aab1 <- aab1[,c(21, 2:20)]
aab1 <- t(aab1)
colnames(aab1) <- aab1[1,]
aab1 <- aab1[-1,]
#make numeric
aab1 <- as.data.frame(aab1)
aab1[,1:473] <- lapply(aab1[,1:473], as.numeric)
aab1 <- as.data.frame(aab1)

a1 <- as.data.frame(colnames(rna2))
names(a1) <- 'id'
a2 <- as.data.frame(colnames(prot1))
names(a2) <- 'id'
a3 <- as.data.frame(colnames(metab1))
names(a3) <- 'id'
a4 <- as.data.frame(colnames(c2))
names(a4) <- 'id'
a5 <- as.data.frame(colnames(aab1))
names(a5) <- 'id'


cols <- as.data.frame(unique(rbind(a1,a2,a3,a4,a5)))
cols <- cols %>% filter(id != "")

#now add the rows with missigness
add_missing_columns <- function(df, ids) {
  for (id in ids) {
    if (!(id %in% colnames(df))) {
      df[[id]] <- NA  
    }
  }
  return(df)
}

# Get the IDs from `cols`
ids_to_add <- cols$id

# Apply the function to each of the dataframes
rna2 <- add_missing_columns(rna2, ids_to_add)
prot1 <- add_missing_columns(prot1, ids_to_add)
metab1 <- add_missing_columns(metab1, ids_to_add)
ab1 <- add_missing_columns(aab1, ids_to_add)
ciber1 <- add_missing_columns(c2, ids_to_add)

####################################################################
#now before filtering, only select relevant metabolites & proteins!
metabolites_conversion_file <- readRDS("/mnt/session_data/mounted-data-readonly/metabolites_conversion_file.rds")
proteomics_proteins_to_keep <- readRDS("/mnt/session_data/mounted-data-readonly/proteomics_proteins_to_keep.rds")

prot1$prots <- rownames(prot1)
#filter proteins
prot1 <- prot1 %>% filter(prots %in% proteomics_proteins_to_keep$id_protein)
#this substantially reduces the data, but corresponds! 
prot1 <- prot1 %>% select(-prots)

#also do this for the metabolites
#renaming has luckily been done now! 
metab1$metabs <- rownames(metab1)

#PROBLEM! Some vars were renamed!!
metab1 <- metab1 %>% filter(metabs %in% metabolites_conversion_file$id_metabolite)
missing_ids <- as.data.frame(setdiff(metabolites_conversion_file$id_metabolites, metab1$metabs))

metab_conv <- unique(metabolites_conversion_file$Metabolite)
#no missing, so that makes sense! 

#so perhaps filter here on missing in other file
#From paper Rhodes et al. 2017. Circ: 625 total participants
#so should be able to calculate the proportion of missing ones! 
metabolites_conversion_file$missing_percent <- (as.numeric(metabolites_conversion_file$Missing)/625)*100
metab_not_keep <- metabolites_conversion_file %>% filter(missing_percent >5)
metab1 <- metab1 %>% filter(!metabs %in% metab_not_keep$id_metabolite)

#for now proceed, check with RJD later
metab1 <- metab1 %>% select(-metabs)


###################
#CONTINUE HERE
##################

#column names need to be in the same order to work....
desired_order1 <- c(colnames(metab1))

rna2 <- rna2 %>% select(all_of(desired_order1))
prot2 <- prot1 %>% select(all_of(desired_order1))
ciber2 <- ciber1 %>% select(all_of(desired_order1))
ab2 <- ab1 %>% select(all_of(desired_order1))



#the data is in shape, now filter on variance!
add_variance_column <- function(df) {
  variance_per_row <- apply(df, 1, function(x) var(x, na.rm = TRUE))
  df$Variance <- variance_per_row
  return(df)
}

rna2 <- add_variance_column(rna2)
prot1 <- add_variance_column(prot1)
metab1 <- add_variance_column(metab1)

#EdB stopped here on 02/10/24
#EdB continued here on 5/11/24 --> used previous filters
#EdB update 30/01/25: now with new filters though 
rna2 <- rna2 %>% arrange(desc(Variance)) %>% head(4524)
prot2 <- prot1 %>% arrange(desc(Variance)) %>% head(2174)
metab1$metabolite_org <- rownames(metab1)

metab2 <- metab1 %>% arrange(desc(Variance)) %>% head(329)

rna2 <- rna2 %>% select(!Variance)
prot2 <- prot2 %>% select(!Variance)
metab2 <- metab2 %>% select(!Variance)

rna2 <- as.matrix(rna2)
prot2 <- as.matrix(prot2)
metab2 <- as.matrix(metab2)
ab2 <- as.matrix(ab2)
ciber2 <- as.matrix(ciber2)

#need to reorder
desired_order1 <- c(colnames(metab2))
rna2 <- as.data.frame(rna2)
prot2 <- as.data.frame(prot2)
ab2 <- as.data.frame(ab2)
ciber2 <- as.data.frame(ciber2)
rna2 <- rna2 %>% select(all_of(desired_order1))
prot2 <- prot2 %>% select(all_of(desired_order1))
ciber2 <- ciber2 %>% select(all_of(desired_order1))
ab2 <- ab2 %>% select(all_of(desired_order1))

rna2 <- as.matrix(rna2)
prot2 <- as.matrix(prot2)
ab2 <- as.matrix(ab2)
ciber2 <- as.matrix(ciber2)
#now make a list
mofa_matrix_full <- list(rna2, prot2, metab2, ab2, ciber2)

setwd('~/session_data')
save(mofa_matrix_full, file = "MOFA_matrix_complete_including_missigness_after_filtering_30Jan25.RData")


#################################################################################################################
# RJ checked w/c 29.9.26. All notes and comments start # RJ
#Defintive MOFA - model 1
#this is CIBER adjusted including leukocytes with no missingness! 

#EdB started here 09/12/24
library(tidyverse)
library(survival)
library(survminer)
library(MOFA2)
library(fgsea)
library(msigdbr)
library(data.table)
library(pheatmap)
library(Publish)
getwd()
library(dplyr)
########################################################################################
#get mortality data ready
clin_dat_at_diagnosis_clean <- readRDS("/mnt/session_data/mounted-data-readonly/PAH-OpenClinica/20240726-freeze/clin_dat_at_diagnosis_clean.rds")
follow_up_data <- readRDS("/mnt/mounted-data/PAH-OpenClinica/20240726-freeze/follow_up_data.rds")

time_df<- clin_dat_at_diagnosis_clean %>% select(id, diagnosis_date)
fudf <- follow_up_data %>% filter(visit == 1) %>% select(id, visit_date)
time_df <- merge(time_df, fudf, by='id')
time_df$visit_date <- as.Date(time_df$visit_date)
time_df$time_to_sampling <- as.numeric(time_df$visit_date - time_df$diagnosis_date)/365.24
time_df <- time_df[,c(1,3,4)]

mort_df <-  clin_dat_at_diagnosis_clean %>% dplyr::select('id', 'diagnosis_verified', 'sex', 'age_diagnosis', 'bs_bmi', 'diagnosis_date', 'sub_cause', 'sub_date', 'hb_pvr_calc', 'hb_rap_m', 'egfr_mdrd', 'cbt_inflammation_crp_mgpl')
mort_df <- merge(mort_df, time_df, by= 'id')
mort_df$sub_date <- as.Date(mort_df$sub_date)
mort_df$diagnosis_date <- as.Date(mort_df$diagnosis_date)
mort_df <- unique(mort_df)
Cenus_1Dec_2023_UK_PAH_cohort <- readRDS("/mnt/mounted-data/Cenus_1Dec_2023_UK_PAH_cohort.rds")

mort_df <- merge(mort_df, Cenus_1Dec_2023_UK_PAH_cohort, by='id')

mort_df$census_cause <- ifelse(mort_df$census_cause == 'DEATH', 'death', mort_df$census_cause)

# RJ contains children and mixed diagnoses at this point

#make sure sub_date is complete and format the data
mort_df$sub_date <- as.Date(mort_df$death_or_census_date)
mort_df$sub_date2 <- fifelse(is.na(mort_df$sub_date), ymd("2023-12-01"), mort_df$sub_date)
mort_df$surv_time <- (mort_df$sub_date2 - mort_df$diagnosis_date)
mort_df[c('surv_time', 'day')] <- str_split_fixed(mort_df$surv_time, ' ', 2)
mort_df$surv_time <- as.numeric(mort_df$surv_time)
mort_df$surv_time <- mort_df$surv_time/365.2422
mort_df$event <- ifelse(mort_df$census_cause == 'death', 1, 0)
mort_df$event <- ifelse(is.na(mort_df$event), 0, mort_df$event)
mort_df$diagnosis_verified <- as.factor(mort_df$diagnosis_verified)
mort_df$sex <- as.factor(mort_df$sex)
mort_df <- mort_df %>% filter(!is.na(surv_time))
#remove negative survival times (due to census date)
mort_df <- mort_df %>% filter(surv_time >=0)
mort_df <- unique(mort_df)

#censor at 15 years for cohort
mort_df$event <- ifelse(mort_df$surv_time >15 & mort_df$event == 1, 0, mort_df$event)
mort_df$surv_time <- ifelse(mort_df$surv_time >15, 15, mort_df$surv_time)

mort_df_for_later <- mort_df

clin <- clin_dat_at_diagnosis_clean %>% select(id, age_diagnosis, sex, bs_bmi, diagnosis_verified, functional_class, cbt_card_ntprobnp_ngpl, cbt_inflammation_crp_mgpl, cbt_thyr_freet4_pmolpl, cbt_thyr_tsh_mupl, cbt_haem_platelets_x10e9pl, cbt_haem_redcell_dist_pct, cbt_haem_wbc_x10e9pl, cbt_haem_hb_gpl, hb_pvr_calc, hb_rap_m, hb_pap_m, hb_cardiac_output_value_1, hb_cardiac_index_value_1, ep_1_distance_meters, hb_pawp_m, lf_fev1_pc, lf_fvc_pc, lf_kco_pc)
######################################################################
#load in the data
#first full matrix

#NOTE: this is variability based filtered! 
load("/mnt/session_data/MOFA_matrix_complete_including_missigness_after_filtering_30Jan25.RData")
#omit this load below, as this is the old dataframe
#load("/mnt/session_data/MOFA_matrix_complete_including_missigness.RData")
load("/mnt/session_data/MOFA_matrix_list_filtered_pts.RData")
Filtered_for_MOFA_TPM_res_PC_WBC_HV_PTs <- readRDS("/mnt/session_data/Filtered_for_MOFA_TPM_res_PC_WBC_HV_PTs.rds")
DOUBLE_check_patient_WBC_PC_residuals <- readRDS("/mnt/session_data/DOUBLE_check_patient_WBC_PC_residuals.rds")
######################################################
#NOTE HERE: THIS IS RE-DONE WITH THE MOFA!
#POTENTIALLY CHANGE THIS BACK (where transcript 
#is selected and centering + scaling applied
#also now turned off scaling the views in the MOFA too
######################################################

a1 <- as.data.frame(filtered_matrix_mofa[[1]])

#transcript <- as.data.frame(mofa_matrix_full[[1]])
transcript <- DOUBLE_check_patient_WBC_PC_residuals
#transcript <- Filtered_for_MOFA_TPM_res_PC_WBC_HV_PTs
#write_rds(transcript, file="residualised_for_leukocytes_andPCs_RNA_seq.rds")


proteome <- as.data.frame(mofa_matrix_full[[2]])
metab <- as.data.frame(mofa_matrix_full[[3]])
aab <- as.data.frame(mofa_matrix_full[[4]])
ciber <- as.data.frame(mofa_matrix_full[[5]])

# RJ fixing nomenclature
# AABs change from _norm to _AAB
rownames(aab) <- sub("_norm$", "_AAB", rownames(aab))

rownames(aab)[rownames(aab) == "Centromere.Protein.B..CENP.B._AAB"] <- "Centromere.Protein.B_AAB"
rownames(aab)[rownames(aab) == "Pyruvate.Dehydrogenase.from.porcine.hear_AAB"] <- "Pyruvate.Dehydrogenase_AAB"
rownames(aab)[rownames(aab) == "Scl.34..Fibrillarin.his.tag._AAB"] <- "Scl.34_AAB"
rownames(aab)[rownames(aab) == "Smith.antigen..Sm._AAB"] <- "Smith.antigen_AAB"

#make leucoctes nicer
rownames(ciber)[rownames(ciber) ==  "ciber.B.cells.naive"] <- "B.cells_Naive"
rownames(ciber)[rownames(ciber) == "ciber.B.cells.memory"] <- "B.cells_Memory"
rownames(ciber)[rownames(ciber) == "ciber.Plasma.cells"] <- "Plasma cells"
rownames(ciber)[rownames(ciber) == "ciber.T.cells.CD8"] <- "CD8+T.cells"
rownames(ciber)[rownames(ciber) == "ciber.T.cells.CD4.naive"] <- "CD4+T.cells_Naive"
rownames(ciber)[rownames(ciber) == "ciber.T.cells.CD4.memory.resting"] <- "Memory_CD4+T.cells_resting"
rownames(ciber)[rownames(ciber) == "ciber.T.cells.CD4.memory.activated"] <- "Memory_CD4+T.cells_activated"
rownames(ciber)[rownames(ciber) == "ciber.T.cells.gamma.delta"] <- "Gamma_delta_T.cells"
rownames(ciber)[rownames(ciber) == "ciber.NK.cells.resting"] <- "NK.cells_resting"
rownames(ciber)[rownames(ciber) == "ciber.NK.cells.activated"] <- "NK.cells_activated"
rownames(ciber)[rownames(ciber) == "ciber.Monocytes"] <- "Monocytes"
rownames(ciber)[rownames(ciber) == "ciber.Macrophages.M0"] <- "M0_Macrophages"
rownames(ciber)[rownames(ciber) == "ciber.Macrophages.M1"] <- "M1_Macrophages"
rownames(ciber)[rownames(ciber) == "ciber.Dendritic.cells.resting"] <- "Dendritic.cells_resting"
rownames(ciber)[rownames(ciber) == "ciber.Dendritic.cells.activated"] <- "Dendritic.cells_activated"
rownames(ciber)[rownames(ciber) == "ciber.Mast.cells.resting"] <- "Mast.cells_resting"
rownames(ciber)[rownames(ciber) == "ciber.Eosinophils"] <- "Eosinophils"
rownames(ciber)[rownames(ciber) == "ciber.Neutrophils"] <- "Neutrophils"


t1 <- transcript %>% select(all_of(names(a1)))
t1 <- as.data.frame(t(t1))
t1[,1:4524] <- lapply(t1[,1:4524], as.numeric)
t1 <- scale(t1, scale=TRUE, center=TRUE)
t1 <- as.data.frame(t(t1))
t1 <- as.matrix(t1)

p1 <- proteome %>% select(all_of(names(a1)))
p1 <- as.data.frame(t(p1))
p1 <- scale(p1, scale=TRUE, center=TRUE)
p1 <- as.data.frame(t(p1))
p1 <- as.matrix(p1)

m1 <- metab %>% select(all_of(names(a1)))
m1 <- as.data.frame(t(m1))
m1 <- scale(m1, scale=TRUE, center=TRUE)
m1 <- as.data.frame(t(m1))
m1 <- as.matrix(m1)

aab1 <- aab %>% select(all_of(names(a1)))
aab1 <- as.data.frame(t(aab1))
aab1 <- scale(aab1, scale=TRUE, center=TRUE)
aab1 <- as.data.frame(t(aab1))
aab1 <- as.matrix(aab1)
aab1 <- as.matrix(aab1)

c1 <- ciber %>% select(all_of(names(a1)))
c1 <- as.data.frame(t(c1))
c1 <- scale(c1, scale=TRUE, center=TRUE)
c1 <- as.data.frame(t(c1))
c1 <- as.matrix(c1)

#this has recuded the number to 179 for each
#make it into a list
corrected_omics_list <- list(t1, p1, m1, aab1, c1)


##################################################################################################
#do the MOFA

#get the code working for the initial MOFA
MOFAobject <- create_mofa(corrected_omics_list)

#gives warning of duplicates features (probably proteomics and RNA-seq)
#should be okay as MOFA renames these
#NOTE: consider log-transforming/normalising the CIBER (albeit MOFA will scale too)
#CIBER was log-transformed for correction for residuals

#visualise MOFA
plot_data_overview(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)
head(model_opts)
#change below for different number of factors (normally 15)
#model_opts$num_factors <- 10
train_opts <- get_default_training_options(MOFAobject)
head(train_opts)
#default seed is 42 (no need to change for now)

data_opts <- get_default_data_options(MOFAobject)
head(data_opts)
#I do want to scale the views as everything is on a differing scale!
###########################################################################
# turned off scaling here now, as already scald + centred manually before #
###########################################################################
#data_opts$scale_views <- TRUE

#now train the MOFA model!!
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

outfile = file.path(getwd(),"model.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)

#Warning message:
#  In .quality_control(object, verbose = verbose) :
#  Factor(s) 2 are strongly correlated with the total number of expressed features for at least one of your omics. Such factors appear when there are differences in the total 'levels' between your samples, *sometimes* because of poor normalisation in the preprocessing steps.
#this error has occured before

#now process the MOFA
#get overview of data
plot_data_overview(MOFAobject.trained)

#also visualise differently

df_vis <- data.frame(Modality = c('Transcriptomics', 'Metabolomics', 'Proteomics', 'Auto-antibodies', 'estimated leukocytes (CIBERSORT)'),
                     vals = c(4524, 329, 2174, 19, 18))



df_vis$Modality <- factor(df_vis$Modality, 
                          levels = c("Transcriptomics", "Proteomics", "Metabolomics", "Auto-antibodies", "estimated leukocytes (CIBERSORT)"))

# Define custom colors
custom_colors <- c("Transcriptomics" = "salmon2", 
                   "CIBERSORT estimated leukocytes" = "grey", 
                   "Proteomics" = "darkgoldenrod1",  
                   "Metabolomics" = "blue", 
                   "Auto-antibodies" = "darkgreen")



plot_vars <- ggplot(df_vis, aes(x=Modality, y=vals, fill=Modality)) +
  geom_bar(stat = "identity") +
  ggtitle('Barchart of number of variables per modality') +
  geom_text(aes(label = vals), vjust = -0.5, size = 5) +
  xlab('Modality') +
  scale_fill_manual(values = custom_colors) +
  ylab('Number of variables') +
  scale_y_continuous(limits = c(0, 5000)) +
  theme_pubr()

print(plot_vars)

plot_vars <- ggplot(df_vis, aes(x = Modality, y = vals, fill = Modality)) +
  geom_bar(stat = "identity") +
  ggtitle('Barchart of number of variables per modality') +
  geom_text(aes(label = vals, 
                vjust = ifelse(vals > 200, 1.5, -0.5),
                color = ifelse(vals > 200, "white", "black")),
            size = 5, fontface = "bold") +
  scale_color_identity() +
  xlab('Modality') +
  scale_fill_manual(values = custom_colors) +
  ylab('Number of variables') +
  theme_pubr() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


#add metadata
MOFAobject.trained@samples_metadata
sum(MOFAobject.trained@dimensions$N)

#also get the factor scores
factor_scores <- get_factors(MOFAobject.trained)
factor_scores_definitive <- as.data.frame(factor_scores)
factor_scores_definitive$samples <- rownames(factor_scores_definitive)

#write_rds(factor_scores_definitive, file='definitive_mofa_model_May25_PC_WBC_resids.rds')

views_names(MOFAobject.trained) <- c("Transcriptomics", "Proteomics", "Metabolomics", "Autoantibodies", "Leukocytes")
varp <- plot_variance_explained(MOFAobject.trained, x="view", y="factor") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Generate the plot
varp <- plot_variance_explained(MOFAobject.trained, x = "view", y = "factor")
print(varp)

r2 <- get_variance_explained(MOFAobject.trained)

# 1) Total variance explained per view (what the top bar plot shows),
#    i.e. the sum across all factors for each modality
r2$r2_total
# list of one element per group (here "group1"), each a named vector

# 2) Variance explained per factor per view (what the heatmap shows)
r2$r2_per_factor


##########################################################################
#Eck addition, plot with varaibles to take forward

# Get the data (run after renaming the views)
r2_factor_mat <- get_variance_explained(MOFAobject.trained)$r2_per_factor[[1]]

# Long format, with a flag for >1%
r2_df <- as.data.frame(r2_factor_mat) %>%
  rownames_to_column("Factor") %>%
  pivot_longer(-Factor, names_to = "View", values_to = "R2") %>%
  mutate(
    Factor = factor(Factor, levels = rev(rownames(r2_factor_mat))),
    View   = factor(View, levels = colnames(r2_factor_mat)),
    flag   = R2 > 1,
    # numeric positions of each tile, needed for the segments
    x = as.numeric(View),
    y = as.numeric(Factor)
  )

r2_df <- r2_df %>% mutate(across(where(is.numeric), ~ ifelse(.x < 1, 0, .x)))

# Subset to flagged tiles only
cross_df <- dplyr::filter(r2_df, flag)

varp_new <- ggplot(r2_df, aes(x = View, y = Factor, fill = R2)) +
  geom_tile(colour = "grey85") +
  scale_fill_gradientn(
    colours = c("white", "#FFF3B0", "#FDB863", "#E66101", "#B2182B", "#67001F"),
    values  = scales::rescale(c(0, 2, 5, 10, 15, 20)),
    limits  = c(0, 20),
    oob     = scales::squish,
    name    = "Var. (%)"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid  = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)
  )

print(varp_new)
#########################################################################################


# Extract actual factor levels from the plot
factors_used <- ggplot_build(varp)$layout$panel_params[[1]]$y$get_labels()

# Reverse their order
varp <- varp +
  scale_y_discrete(limits = rev(factors_used)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



# Display the plot
varp


varp <- varp +
  scale_fill_gradientn(
    colours = c(
      "white",    # 0%
      "#FFF3B0",  # 2%
      "#FDB863",  # 5%
      "#E66101",  # 10%
      "#B2182B",  # 15%
      "#67001F"   # 20%
    ),
    values = scales::rescale(c(0, 2, 5, 10, 15, 20)),
    limits = c(0, 20),
    name = "Var. (%)"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 1,
      vjust = 0.5
    )
  )

print(varp)

varp <- varp +
  scale_y_discrete(limits = rev(factors_used)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


########################################################################
#attempt combined plot:
df_vis$Modality <- recode(df_vis$Modality,
                          "Auto-antibodies"                      = "Autoantibodies",
                          "estimated leukocytes (CIBERSORT)"     = "Leukocytes"
)


modality_order <- c("Transcriptomics", "Proteomics", "Metabolomics", 
                    "Autoantibodies", "Leukocytes")

# Rj - didnt bother with this bit

df_vis$Modality <- factor(df_vis$Modality, levels = modality_order)

#Update custom_colors names to match too
names(custom_colors) <- recode(names(custom_colors),
                               "Auto-antibodies"                      = "Autoantibodies",
                               "estimated leukocytes (CIBERSORT)"     = "Leukocytes"
)


plot_vars <- ggplot(df_vis, aes(x = Modality, y = vals, fill = Modality)) +
  geom_bar(stat = "identity") +
  ggtitle('Barchart of number of variables per modality') +
  geom_text(aes(label = vals, 
                vjust = ifelse(vals > 330, 1.5, -0.5),
                color = ifelse(vals > 330, "white", "black")),
            size = 5, fontface = "bold") +
  scale_color_identity() +
  xlab('Modality') +
  scale_fill_manual(values = custom_colors) +
  ylab('Number of variables') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +  # extra headroom on top
  theme_pubr() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5)  # trim bottom margin
  )


varp <- varp +
  scale_y_discrete(limits = rev(factors_used)) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 13),
    axis.text.y = element_text(size = 13),
    axis.title = element_text(size = 14),
    plot.margin = margin(t = 0, r = 5, b = 5, l = 5)
  )


library(cowplot)
combined_plot <- plot_grid(
  plot_vars, varp,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1.6))

combined_plot
#########################################################################
#also get weights
weights <- get_weights(MOFAobject.trained, views = "all", factors = "all")
rna_weight <- as.data.frame(weights$Transcriptomics)
proteome_weight <- as.data.frame(weights$Proteomics)
metab_weight <- as.data.frame(weights$Metabolomics)
aab_weight <- as.data.frame(weights$Autoantibodies)
leuko_weight <- as.data.frame(weights$Leukocytes)

#make CSV files for export
write.csv2(rna_weight, file='MOFA_trascriptomic_weights_6Nov_25.csv')
write.csv2(proteome_weight, file='MOFA_proteome_weights_6Nov_25.csv')
write.csv2(metab_weight, file='MOFA_metabolomic_weights_6Nov_25.csv')
write.csv2(aab_weight, file='MOFA_autoantibodies_weights_6Nov_25.csv')
write.csv2(leuko_weight, file='MOFA_leukocytes_weights_6Nov_25.csv')

lapply(proteome_weight, hist)
lapply(metab_weight, hist)



#now investigate factors
plot_top_weights(MOFAobject.trained,
                 view = "Transcriptomics",
                 factor = 1,
                 nfeatures = 25
)


plot_top_weights(MOFAobject.trained,
                 view = "Autoantibodies",
                 factor = 1,
                 nfeatures = 25
)


plot_top_weights(MOFAobject.trained,
                 view = "Proteomics",
                 factor = 1,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Proteomics",
                 factor = 2,
                 nfeatures = 25
)
plot_top_weights(MOFAobject.trained,
                 view = "Metabolomics",
                 factor = 2,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Leukocytes",
                 factor = 2,
                 nfeatures = 25
)


plot_top_weights(MOFAobject.trained,
                 view = "Proteomics",
                 factor = 3,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Metabolomics",
                 factor = 3,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Leukocytes",
                 factor = 3,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Metabolomics",
                 factor = 3,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Proteomics",
                 factor = 4,
                 nfeatures = 25
)


plot_top_weights(MOFAobject.trained,
                 view = "Proteomics",
                 factor = 5,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Leukocytes",
                 factor = 5,
                 nfeatures = 25)

plot_top_weights(MOFAobject.trained,
                 view = "Proteomics",
                 factor = 6,
                 nfeatures = 25
)


plot_top_weights(MOFAobject.trained,
                 view = "Metabolomics",
                 factor = 6,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Leukocytes",
                 factor = 6,
                 nfeatures = 25)



plot_top_weights(MOFAobject.trained,
                 view = "Transcriptomics",
                 factor = 7,
                 nfeatures = 25
)
plot_top_weights(MOFAobject.trained,
                 view = "Proteomics",
                 factor = 7,
                 nfeatures = 25
)
plot_top_weights(MOFAobject.trained,
                 view = "Metabolomics",
                 factor = 7,
                 nfeatures = 25
)
plot_top_weights(MOFAobject.trained,
                 view = "Autoantibodies",
                 factor = 7,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Transcriptomics",
                 factor = 8,
                 nfeatures = 25
)


plot_top_weights(MOFAobject.trained,
                 view = "Transcriptomics",
                 factor = 9,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Proteomics",
                 factor = 10,
                 nfeatures = 25
)
plot_top_weights(MOFAobject.trained,
                 view = "Metabolomics",
                 factor = 10,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Transcriptomics",
                 factor = 11,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Transcriptomics",
                 factor = 12,
                 nfeatures = 25
)


plot_top_weights(MOFAobject.trained,
                 view = "Transcriptomics",
                 factor = 13,
                 nfeatures = 25
)


plot_top_weights(MOFAobject.trained,
                 view = "Transcriptomics",
                 factor = 14,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Transcriptomics",
                 factor = 15,
                 nfeatures = 25
)

#visualise factor 13 weights
fact13 <- rna_weight %>% filter(abs(Factor13) >0.5)
rownames(fact13)[rownames(fact13) == "OAS1_view_1"] <- "OAS1"
rownames(fact13)[rownames(fact13) == "DDX58_view_1"] <- "DDX58"
fact13$gene <- rownames(fact13)


library(ggplot2)
library(ggpubr)


ggplot(fact13, aes(x = reorder(gene, -Factor13), y = Factor13)) +
  geom_bar(stat = 'identity', fill = "#C8A2C8") +
  xlab('Gene') +
  ylab('Weighting of gene') +
  ggtitle('Gene weight for each transcript’s contribution to factor 13') +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


######################################################################################
#quickly check survival!
#get survival differences
#merge with MOFA datamort_df#merge with MOFA data
mdf <- merge(mort_df, factor_scores_definitive, by.x='id', by.y='samples')
hist(mdf$time_to_sampling, breaks=50)
#as censored at 15 years do 15 year max time to sampling
#mdf$time_to_sampling <- ifelse(mdf$time_to_sampling >15, 15, mdf$time_to_sampling)
hist(log(mdf$hb_rap_m +10))
hist(log(mdf$bs_bmi))
lm1 <- lm(mdf$group1.Factor13 ~ mdf$age_diagnosis + mdf$sex + mdf$time_to_sampling + mdf$hb_rap_m + mdf$hb_pvr_calc + mdf$bs_bmi +mdf$cbt_inflammation_crp_mgpl)
summary(lm1)

lm2 <- lm(mdf$time_to_sampling ~ mdf$age_diagnosis + mdf$sex + mdf$hb_rap_m + mdf$hb_pvr_calc + mdf$bs_bmi +mdf$cbt_inflammation_crp_mgpl)
lm3 <- lm(mdf$time_to_sampling ~ mdf$age_diagnosis + mdf$sex + mdf$hb_rap_m + mdf$hb_pvr_calc + mdf$bs_bmi +mdf$cbt_inflammation_crp_mgpl +mdf$group1.Factor13)

summary(lm2)
summary(lm3)

sobj <- Surv(mdf$surv_time, mdf$event, type='right')
cox <- coxph(sobj ~ group1.Factor1 + group1.Factor2 + group1.Factor3 + group1.Factor4 + group1.Factor5 + group1.Factor6 + group1.Factor7 + group1.Factor8 + group1.Factor9 + group1.Factor10 + group1.Factor11 + group1.Factor12 + group1.Factor13+ group1.Factor14 + group1.Factor15, data=mdf)
ggforest(cox, data=mdf)

#RJ - at this point the cox does not match the one in the paper

#extra check without factors
cox <- coxph(sobj ~ age_diagnosis + sex + hb_pvr_calc + hb_rap_m + bs_bmi + cbt_inflammation_crp_mgpl + time_to_sampling, data=mdf)
ggforest(cox, data=mdf)

#now select those of interest
cox <- coxph(sobj ~ group1.Factor1 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + bs_bmi + cbt_inflammation_crp_mgpl, data=mdf)
#cox <- coxph(sobj ~ group1.Factor1 + age_diagnosis + sex + time_to_sampling, data=mdf)
ggforest(cox, data=mdf)
summary_cox <- summary(cox)

#Extract coefs
hr <- summary_cox$coefficients[, "exp(coef)"]
ci_lower <- summary_cox$conf.int[, "lower .95"]
ci_upper <- summary_cox$conf.int[, "upper .95"]
pval <- summary_cox$coefficients[, "Pr(>|z|)"]

#get in df
cox_f1 <- data.frame(
  Variable = rownames(summary_cox$coefficients),
  HR = round(hr, 3),
  CI_lower_95 = round(ci_lower, 3),
  CI_upper_95 = round(ci_upper, 3),
  p_value = signif(pval, 3))



cox <- coxph(sobj ~ group1.Factor2 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + bs_bmi + cbt_inflammation_crp_mgpl, data=mdf)
#cox <- coxph(sobj ~ group1.Factor2 + age_diagnosis + sex + time_to_sampling, data=mdf)
ggforest(cox, data=mdf)
summary_cox <- summary(cox)

#Extract coefs
hr <- summary_cox$coefficients[, "exp(coef)"]
ci_lower <- summary_cox$conf.int[, "lower .95"]
ci_upper <- summary_cox$conf.int[, "upper .95"]
pval <- summary_cox$coefficients[, "Pr(>|z|)"]

#get in df
cox_f2 <- data.frame(
  Variable = rownames(summary_cox$coefficients),
  HR = round(hr, 3),
  CI_lower_95 = round(ci_lower, 3),
  CI_upper_95 = round(ci_upper, 3),
  p_value = signif(pval, 3))

cox <- coxph(sobj ~  group1.Factor3 + age_diagnosis + sex+ hb_pvr_calc + hb_rap_m + bs_bmi + cbt_inflammation_crp_mgpl, data=mdf)
#cox <- coxph(sobj ~ group1.Factor3 + age_diagnosis + sex + time_to_sampling, data=mdf)
ggforest(cox, data=mdf)
summary_cox <- summary(cox)

#Extract coefs
hr <- summary_cox$coefficients[, "exp(coef)"]
ci_lower <- summary_cox$conf.int[, "lower .95"]
ci_upper <- summary_cox$conf.int[, "upper .95"]
pval <- summary_cox$coefficients[, "Pr(>|z|)"]

#get in df
cox_f3 <- data.frame(
  Variable = rownames(summary_cox$coefficients),
  HR = round(hr, 3),
  CI_lower_95 = round(ci_lower, 3),
  CI_upper_95 = round(ci_upper, 3),
  p_value = signif(pval, 3))

cox <- coxph(sobj ~ group1.Factor4 + age_diagnosis + sex+ hb_pvr_calc + hb_rap_m + bs_bmi + cbt_inflammation_crp_mgpl, data=mdf)
#cox <- coxph(sobj ~ group1.Factor4 + age_diagnosis + sex + time_to_sampling, data=mdf)
ggforest(cox, data=mdf)
summary_cox <- summary(cox)

#Extract coefs
hr <- summary_cox$coefficients[, "exp(coef)"]
ci_lower <- summary_cox$conf.int[, "lower .95"]
ci_upper <- summary_cox$conf.int[, "upper .95"]
pval <- summary_cox$coefficients[, "Pr(>|z|)"]

#get in df
cox_f4 <- data.frame(
  Variable = rownames(summary_cox$coefficients),
  HR = round(hr, 3),
  CI_lower_95 = round(ci_lower, 3),
  CI_upper_95 = round(ci_upper, 3),
  p_value = signif(pval, 3))

cox <- coxph(sobj ~ group1.Factor5 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + bs_bmi + cbt_inflammation_crp_mgpl, data=mdf)
#cox <- coxph(sobj ~ group1.Factor5 + age_diagnosis + sex + time_to_sampling, data=mdf)
ggforest(cox, data=mdf)
summary_cox <- summary(cox)

#Extract coefs
hr <- summary_cox$coefficients[, "exp(coef)"]
ci_lower <- summary_cox$conf.int[, "lower .95"]
ci_upper <- summary_cox$conf.int[, "upper .95"]
pval <- summary_cox$coefficients[, "Pr(>|z|)"]

#get in df
cox_f5 <- data.frame(
  Variable = rownames(summary_cox$coefficients),
  HR = round(hr, 3),
  CI_lower_95 = round(ci_lower, 3),
  CI_upper_95 = round(ci_upper, 3),
  p_value = signif(pval, 3))

cox <- coxph(sobj ~ group1.Factor6 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + bs_bmi + cbt_inflammation_crp_mgpl, data=mdf)
#cox <- coxph(sobj ~ group1.Factor6 + age_diagnosis + sex + time_to_sampling, data=mdf)
ggforest(cox, data=mdf)
summary_cox <- summary(cox)

#Extract coefs
hr <- summary_cox$coefficients[, "exp(coef)"]
ci_lower <- summary_cox$conf.int[, "lower .95"]
ci_upper <- summary_cox$conf.int[, "upper .95"]
pval <- summary_cox$coefficients[, "Pr(>|z|)"]

#get in df
cox_f6 <- data.frame(
  Variable = rownames(summary_cox$coefficients),
  HR = round(hr, 3),
  CI_lower_95 = round(ci_lower, 3),
  CI_upper_95 = round(ci_upper, 3),
  p_value = signif(pval, 3))

cox <- coxph(sobj ~ group1.Factor7 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m +bs_bmi + cbt_inflammation_crp_mgpl, data=mdf)
#cox <- coxph(sobj ~ group1.Factor7 + age_diagnosis + sex + time_to_sampling, data=mdf)
ggforest(cox, data=mdf)
summary_cox <- summary(cox)

#Extract coefs
hr <- summary_cox$coefficients[, "exp(coef)"]
ci_lower <- summary_cox$conf.int[, "lower .95"]
ci_upper <- summary_cox$conf.int[, "upper .95"]
pval <- summary_cox$coefficients[, "Pr(>|z|)"]

#get in df
cox_f7 <- data.frame(
  Variable = rownames(summary_cox$coefficients),
  HR = round(hr, 3),
  CI_lower_95 = round(ci_lower, 3),
  CI_upper_95 = round(ci_upper, 3),
  p_value = signif(pval, 3))


cox <- coxph(sobj ~ group1.Factor8 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + bs_bmi+ cbt_inflammation_crp_mgpl, data=mdf)
#cox <- coxph(sobj ~ group1.Factor8 + age_diagnosis + sex + time_to_sampling, data=mdf)
ggforest(cox, data=mdf)
summary_cox <- summary(cox)

#Extract coefs
hr <- summary_cox$coefficients[, "exp(coef)"]
ci_lower <- summary_cox$conf.int[, "lower .95"]
ci_upper <- summary_cox$conf.int[, "upper .95"]
pval <- summary_cox$coefficients[, "Pr(>|z|)"]

#get in df
cox_f8 <- data.frame(
  Variable = rownames(summary_cox$coefficients),
  HR = round(hr, 3),
  CI_lower_95 = round(ci_lower, 3),
  CI_upper_95 = round(ci_upper, 3),
  p_value = signif(pval, 3))


cox <- coxph(sobj ~ group1.Factor9 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m +bs_bmi + cbt_inflammation_crp_mgpl, data=mdf)
#cox <- coxph(sobj ~ group1.Factor9 + age_diagnosis + sex + time_to_sampling, data=mdf)
ggforest(cox, data=mdf)
summary_cox <- summary(cox)

#Extract coefs
hr <- summary_cox$coefficients[, "exp(coef)"]
ci_lower <- summary_cox$conf.int[, "lower .95"]
ci_upper <- summary_cox$conf.int[, "upper .95"]
pval <- summary_cox$coefficients[, "Pr(>|z|)"]

#get in df
cox_f9 <- data.frame(
  Variable = rownames(summary_cox$coefficients),
  HR = round(hr, 3),
  CI_lower_95 = round(ci_lower, 3),
  CI_upper_95 = round(ci_upper, 3),
  p_value = signif(pval, 3))


cox <- coxph(sobj ~ group1.Factor10 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + bs_bmi + cbt_inflammation_crp_mgpl, data=mdf)
#cox <- coxph(sobj ~ group1.Factor10 + age_diagnosis + sex + time_to_sampling, data=mdf)
ggforest(cox, data=mdf)
summary_cox <- summary(cox)

#Extract coefs
hr <- summary_cox$coefficients[, "exp(coef)"]
ci_lower <- summary_cox$conf.int[, "lower .95"]
ci_upper <- summary_cox$conf.int[, "upper .95"]
pval <- summary_cox$coefficients[, "Pr(>|z|)"]

#get in df
cox_f10 <- data.frame(
  Variable = rownames(summary_cox$coefficients),
  HR = round(hr, 3),
  CI_lower_95 = round(ci_lower, 3),
  CI_upper_95 = round(ci_upper, 3),
  p_value = signif(pval, 3))


#extra check!
cox <- coxph(sobj ~ group1.Factor10 + group1.Factor12 +  age_diagnosis + sex + hb_pvr_calc + hb_rap_m + bs_bmi + cbt_inflammation_crp_mgpl, data=mdf)
ggforest(cox, data=mdf) 


cox <- coxph(sobj ~group1.Factor11 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + bs_bmi + cbt_inflammation_crp_mgpl, data=mdf)
#cox <- coxph(sobj ~ group1.Factor11 + age_diagnosis + sex + time_to_sampling, data=mdf)
ggforest(cox, data=mdf)
summary_cox <- summary(cox)

#Extract coefs
hr <- summary_cox$coefficients[, "exp(coef)"]
ci_lower <- summary_cox$conf.int[, "lower .95"]
ci_upper <- summary_cox$conf.int[, "upper .95"]
pval <- summary_cox$coefficients[, "Pr(>|z|)"]

#get in df
cox_f11 <- data.frame(
  Variable = rownames(summary_cox$coefficients),
  HR = round(hr, 3),
  CI_lower_95 = round(ci_lower, 3),
  CI_upper_95 = round(ci_upper, 3),
  p_value = signif(pval, 3))


cox <- coxph(sobj ~ group1.Factor12 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + bs_bmi + cbt_inflammation_crp_mgpl, data=mdf)
#cox <- coxph(sobj ~ group1.Factor12 + age_diagnosis + sex + time_to_sampling, data=mdf)
ggforest(cox, data=mdf)
summary_cox <- summary(cox)

#Extract coefs
hr <- summary_cox$coefficients[, "exp(coef)"]
ci_lower <- summary_cox$conf.int[, "lower .95"]
ci_upper <- summary_cox$conf.int[, "upper .95"]
pval <- summary_cox$coefficients[, "Pr(>|z|)"]

#get in df
cox_f12 <- data.frame(
  Variable = rownames(summary_cox$coefficients),
  HR = round(hr, 3),
  CI_lower_95 = round(ci_lower, 3),
  CI_upper_95 = round(ci_upper, 3),
  p_value = signif(pval, 3))

cox <- coxph(sobj ~ group1.Factor13 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + bs_bmi + cbt_inflammation_crp_mgpl, data=mdf)
#cox <- coxph(sobj ~ group1.Factor13 + age_diagnosis + sex + time_to_sampling, data=mdf)
ggforest(cox, data=mdf)
summary(cox)
summary_cox <- summary(cox)

#Extract coefs
hr <- summary_cox$coefficients[, "exp(coef)"]
ci_lower <- summary_cox$conf.int[, "lower .95"]
ci_upper <- summary_cox$conf.int[, "upper .95"]
pval <- summary_cox$coefficients[, "Pr(>|z|)"]

#get in df
cox_f13 <- data.frame(
  Variable = rownames(summary_cox$coefficients),
  HR = round(hr, 3),
  CI_lower_95 = round(ci_lower, 3),
  CI_upper_95 = round(ci_upper, 3),
  p_value = signif(pval, 3))


cox <- coxph(sobj ~ group1.Factor14 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + bs_bmi + cbt_inflammation_crp_mgpl, data=mdf)
#cox <- coxph(sobj ~ group1.Factor14 + age_diagnosis + sex + time_to_sampling, data=mdf)
ggforest(cox, data=mdf)
summary_cox <- summary(cox)

#Extract coefs
hr <- summary_cox$coefficients[, "exp(coef)"]
ci_lower <- summary_cox$conf.int[, "lower .95"]
ci_upper <- summary_cox$conf.int[, "upper .95"]
pval <- summary_cox$coefficients[, "Pr(>|z|)"]

#get in df
cox_f14 <- data.frame(
  Variable = rownames(summary_cox$coefficients),
  HR = round(hr, 3),
  CI_lower_95 = round(ci_lower, 3),
  CI_upper_95 = round(ci_upper, 3),
  p_value = signif(pval, 3))

cox <- coxph(sobj ~ group1.Factor15 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + bs_bmi + cbt_inflammation_crp_mgpl, data=mdf)
#cox <- coxph(sobj ~ group1.Factor15 + age_diagnosis + sex + time_to_sampling, data=mdf)
ggforest(cox, data=mdf)
summary_cox <- summary(cox)

#Extract coefs
hr <- summary_cox$coefficients[, "exp(coef)"]
ci_lower <- summary_cox$conf.int[, "lower .95"]
ci_upper <- summary_cox$conf.int[, "upper .95"]
pval <- summary_cox$coefficients[, "Pr(>|z|)"]

#get in df
cox_f15 <- data.frame(
  Variable = rownames(summary_cox$coefficients),
  HR = round(hr, 3),
  CI_lower_95 = round(ci_lower, 3),
  CI_upper_95 = round(ci_upper, 3),
  p_value = signif(pval, 3))

cox_sum <- rbind(cox_f1[1,], cox_f2[1,],cox_f3[1,],cox_f4[1,],cox_f5[1,],cox_f6[1,],cox_f7[1,],cox_f8[1,],cox_f9[1,],cox_f10[1,],cox_f11[1,],cox_f12[1,],cox_f13[1,],cox_f14[1,],cox_f15[1,])
#Rename vars
cox_sum$Variable <- paste("Factor", seq_len(nrow(cox_sum)))


#plot this
cox_sum$Factor <- factor(cox_sum$Variable, levels = rev(cox_sum$Variable))

#position pvals
x_text_pos <- max(cox_sum$CI_upper_95) * 1.15

# Plot
survp1 <- ggplot(cox_sum, aes(x = HR, y = Factor)) +
  geom_point(size = 3, color = "blue") +  # HR point
  geom_errorbarh(aes(xmin = CI_lower_95, xmax = CI_upper_95), height = 0.2) +  
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +  
  geom_text(aes(x = x_text_pos, label = paste0("p = ", signif(p_value, 3))), 
            size = 3.5, hjust = 0, color = "black") +  
  labs(x = "Hazard Ratio (HR)", y = "Factor", 
       title = "Cox Proportional Hazards Results") +
  theme_pubr() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12)
  ) +
  xlim(min(cox_sum$CI_lower_95) * 0.9, x_text_pos * 1.1)  

survp1


cox_sum$fdr <- p.adjust(cox_sum$p_value, method = 'fdr')

survp1 <- ggplot(cox_sum, aes(x = HR, y = Factor)) +
  geom_point(size = 3, color = "blue") +  # HR point
  geom_errorbarh(aes(xmin = CI_lower_95, xmax = CI_upper_95), height = 0.2) +  
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +  
  geom_text(aes(x = x_text_pos, label = paste0("FDR = ", signif(fdr, 3))), 
            size = 3.5, hjust = 0, color = "black") +  
  labs(x = "Hazard Ratio (HR)", y = "Factor", 
       title = "Cox Proportional Hazards Results") +
  theme_pubr() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12)
  ) +
  xlim(min(cox_sum$CI_lower_95) * 0.9, x_text_pos * 1.1) 

survp1

x_hr  <- max(cox_sum$CI_upper_95) * 1.1
x_fdr <- x_hr * 1.15          

survp2 <- ggplot(cox_sum, aes(x = HR, y = Factor)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbarh(aes(xmin = CI_lower_95, xmax = CI_upper_95), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_text(aes(x = x_hr,
                label = paste0("HR = ", formatC(HR, format = "f", digits = 2))),
            size = 3.5, hjust = 0, color = "black") +
  geom_text(aes(x = x_fdr,
                label = paste0("(FDR = ", formatC(fdr, format = "f", digits = 3), ")"),
                color = fdr < 0.05),
            size = 3.5, hjust = 0) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), guide = "none") +
  labs(x = "Hazard Ratio (HR)",
       title = "Cox Proportional Hazards Results") +
  theme_pubr() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y  = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title   = element_text(size = 14, face = "bold"),
    axis.title   = element_text(size = 12)
  ) +
  xlim(min(cox_sum$CI_lower_95) * 0.9, x_fdr * 1.1)

survp2

write.csv(cox_sum, "RJ_HR for MOFA.csv")
#========================================================
#correlate with original factors
#get factor correlations
definitive_mofa_model_factor_scores <- readRDS("/mnt/session_data/definitive_mofa_model_factor_scores.rds")

# RJ - these definitive MOFA scores are NOT the same as the factor_scores_definitive - is this the one with missingness?
factors_first_mofa <- definitive_mofa_model_factor_scores
#now get correlations between factors
colnames(factors_first_mofa) <- ifelse(colnames(factors_first_mofa) == "samples", "samples", paste0("original_", colnames(factors_first_mofa)))
colnames(factor_scores_definitive) <- ifelse(colnames(factor_scores_definitive) == "samples", "samples", paste0("filtered__", colnames(factor_scores_definitive)))
#now names are changes, plot in a df for correlations
df1 <- merge(factors_first_mofa, factor_scores_definitive, by='samples')
df1 <- df1[,-1]

#now plot!
row_vars <- grep("^original_", colnames(df1), value = TRUE)
col_vars <- setdiff(colnames(df1), row_vars)

# Subset the dataframe
df_rows <- df1[, row_vars]
df_cols <- df1[, col_vars]

# Calculate the correlation matrix between row and column variables
cor_matrix <- cor(df_rows, df_cols, method = "pearson")

# Create the heatmap
pheatmap(cor_matrix, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE)

# RJ - this does not match S4A
#EdB comment: that is correct, plot for S4 is done later downstream 
#=======================================================================================================
#Do RNA-seq pathway analysis
#now do pathway analysis per variable (using fgsea)
library(fgsea)
library(msigdbr)

rna_ranked <- rna_weight
rna_ranked[,1:15] <- lapply(rna_ranked[,1:15], abs)
# Convert the result back to a dataframe (if needed)
rna_ranked_df <- as.data.frame(rna_ranked)
rna_ranked_df$gene <-rownames(rna_ranked_df)

#remove the view_1 so that the gene name is properly recorded
rna_ranked_df[] <- lapply(rna_ranked_df, function(x) if(is.character(x)) sub("_view_1", "", x) else x)
rownames(rna_ranked_df) <- rna_ranked_df$gene

ranked_genes1 <- setNames(rna_ranked_df$Factor1, rna_ranked_df$gene)
ranked_genes2 <- setNames(rna_ranked_df$Factor2, rna_ranked_df$gene)
ranked_genes3 <- setNames(rna_ranked_df$Factor3, rna_ranked_df$gene)
ranked_genes4 <- setNames(rna_ranked_df$Factor4, rna_ranked_df$gene)
ranked_genes5 <- setNames(rna_ranked_df$Factor5, rna_ranked_df$gene)
ranked_genes6 <- setNames(rna_ranked_df$Factor6, rna_ranked_df$gene)
ranked_genes7 <- setNames(rna_ranked_df$Factor7, rna_ranked_df$gene)
ranked_genes8 <- setNames(rna_ranked_df$Factor8, rna_ranked_df$gene)
ranked_genes9 <- setNames(rna_ranked_df$Factor9, rna_ranked_df$gene)
ranked_genes10 <- setNames(rna_ranked_df$Factor10, rna_ranked_df$gene)
ranked_genes11 <- setNames(rna_ranked_df$Factor11, rna_ranked_df$gene)
ranked_genes12 <- setNames(rna_ranked_df$Factor12, rna_ranked_df$gene)
ranked_genes13 <- setNames(rna_ranked_df$Factor13, rna_ranked_df$gene)
ranked_genes14 <- setNames(rna_ranked_df$Factor14, rna_ranked_df$gene)
ranked_genes15 <- setNames(rna_ranked_df$Factor15, rna_ranked_df$gene)

# RJ not sure there was any point in generating all of these if only using a subset 


pathways <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  split(x = .$gene_symbol, f = .$gs_name)


#NB: need to think about minsize for pathways! 
#do factor 1
fgsea_results1 <- fgsea(pathways = pathways, 
                        stats = ranked_genes1, 
                        minSize = 10,  
                        maxSize = 500, scoreType = "pos")


res1 <- fgsea_results1
res1 <- res1 %>% filter(padj <0.05)

ggplot(res1, aes(reorder(pathway, NES), NES)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "FGSEA Results - Factor 1",
    x = "Pathway",
    y = "Normalised Enrichment Score (NES)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),  # Reduce size of y-axis labels (pathways)
    axis.text.x = element_text(size = 8), # Adjust size of x-axis labels if needed
    plot.title = element_text(size = 14, face = "bold")  # Title size and style
  )

fgsea_results7 <- fgsea(pathways = pathways, 
                        stats = ranked_genes7, 
                        minSize = 10,  
                        maxSize = 500, scoreType = "pos")  

res1 <- fgsea_results7
res1 <- res1 %>% filter(padj <0.05)

ggplot(res1, aes(reorder(pathway, NES), NES)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "FGSEA Results - Factor 7",
    x = "Pathway",
    y = "Normalised Enrichment Score (NES)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),  # Reduce size of y-axis labels (pathways)
    axis.text.x = element_text(size = 8), # Adjust size of x-axis labels if needed
    plot.title = element_text(size = 14, face = "bold")  # Title size and style
  )

fgsea_results8 <- fgsea(pathways = pathways, 
                        stats = ranked_genes8, 
                        minSize = 10,  
                        maxSize = 500, scoreType = "pos")  


res1 <- fgsea_results8
res1 <- res1 %>% filter(padj <0.05)

ggplot(res1, aes(reorder(pathway, NES), NES)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "FGSEA Results - Factor 8",
    x = "Pathway",
    y = "Normalised Enrichment Score (NES)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),  # Reduce size of y-axis labels (pathways)
    axis.text.x = element_text(size = 8), # Adjust size of x-axis labels if needed
    plot.title = element_text(size = 14, face = "bold")  # Title size and style
  )

#do factor 9
fgsea_results9 <- fgsea(pathways = pathways, 
                        stats = ranked_genes9, 
                        minSize = 10,  
                        maxSize = 500, scoreType = "pos")  


res1 <- fgsea_results9
res1 <- res1 %>% filter(padj <0.05)

ggplot(res1, aes(reorder(pathway, NES), NES)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "FGSEA Results - Factor 9",
    x = "Pathway",
    y = "Normalised Enrichment Score (NES)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),  # Reduce size of y-axis labels (pathways)
    axis.text.x = element_text(size = 8), # Adjust size of x-axis labels if needed
    plot.title = element_text(size = 14, face = "bold")  # Title size and style
  )


#do factor 11
fgsea_results11 <- fgsea(pathways = pathways, 
                         stats = ranked_genes11, 
                         minSize = 10,  
                         maxSize = 500, scoreType = "pos") 

res1 <- fgsea_results11
res1 <- res1 %>% filter(padj <0.05)

ggplot(res1, aes(reorder(pathway, NES), NES)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "FGSEA Results - Factor 11",
    x = "Pathway",
    y = "Normalised Enrichment Score (NES)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),  # Reduce size of y-axis labels (pathways)
    axis.text.x = element_text(size = 8), # Adjust size of x-axis labels if needed
    plot.title = element_text(size = 14, face = "bold")  # Title size and style
  )

#do factor 12
fgsea_results12 <- fgsea(pathways = pathways, 
                         stats = ranked_genes12, 
                         minSize = 10,  
                         maxSize = 500, scoreType = "pos")  

res1 <- fgsea_results12
res1 <- res1 %>% filter(padj <0.05)

ggplot(res1, aes(reorder(pathway, NES), NES)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "FGSEA Results - Factor 12",
    x = "Pathway",
    y = "Normalised Enrichment Score (NES)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),  # Reduce size of y-axis labels (pathways)
    axis.text.x = element_text(size = 8), # Adjust size of x-axis labels if needed
    plot.title = element_text(size = 14, face = "bold")  # Title size and style
  )

#do factor 13
fgsea_results13 <- fgsea(pathways = pathways, 
                         stats = ranked_genes13, 
                         minSize = 10,  
                         maxSize = 500,
                         scoreType = "pos")  

res1 <- fgsea_results13
res1 <- res1 %>% filter(padj <0.05)

ggplot(res1, aes(reorder(pathway, NES), NES)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "FGSEA Results - Factor 13",
    x = "Pathway",
    y = "Normalised Enrichment Score (NES)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),  # Reduce size of y-axis labels (pathways)
    axis.text.x = element_text(size = 8), # Adjust size of x-axis labels if needed
    plot.title = element_text(size = 14, face = "bold")  # Title size and style
  )

#do factor 14
fgsea_results14 <- fgsea(pathways = pathways, 
                         stats = ranked_genes14, 
                         minSize = 10,  
                         maxSize = 500, scoreType = "pos")  

res1 <- fgsea_results14
res1 <- res1 %>% filter(padj <0.05)

ggplot(res1, aes(reorder(pathway, NES), NES)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "FGSEA Results - Factor 14",
    x = "Pathway",
    y = "Normalised Enrichment Score (NES)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),  # Reduce size of y-axis labels (pathways)
    axis.text.x = element_text(size = 8), # Adjust size of x-axis labels if needed
    plot.title = element_text(size = 14, face = "bold")  # Title size and style
  )

#do factor 15
fgsea_results15 <- fgsea(pathways = pathways, 
                         stats = ranked_genes15, 
                         minSize = 10,  
                         maxSize = 500, scoreType = "pos") 

res1 <- fgsea_results15
res1 <- res1 %>% filter(padj <0.05)

ggplot(res1, aes(reorder(pathway, NES), NES)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "FGSEA Results - Factor 15",
    x = "Pathway",
    y = "Normalised Enrichment Score (NES)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),  # Reduce size of y-axis labels (pathways)
    axis.text.x = element_text(size = 8), # Adjust size of x-axis labels if needed
    plot.title = element_text(size = 14, face = "bold")  # Title size and style
  )

##############################################################################
#Eck edit: ridge plot to improve pathway visualisation --> changed to boxplots for abs()

# ===== Boxplots of FGSEA results =====
box_factors <- c(1, 8, 9, 11, 12, 13, 14, 15)
padj_cut_global <- 0.05

# shared colour scale: NES range across all significant pathways in all factors
all_nes <- unlist(lapply(box_factors, function(k) {
  res <- get(paste0("fgsea_results", k))
  res$NES[!is.na(res$padj) & res$padj < padj_cut_global]
}))
nes_limits <- range(all_nes, na.rm = TRUE)   # c(min, max)
# to anchor the scale at zero instead: nes_limits <- c(0, max(all_nes, na.rm = TRUE))

make_box_plot <- function(k, padj_cut = padj_cut_global) {
  
  res   <- get(paste0("fgsea_results", k))
  stats <- get(paste0("ranked_genes", k))
  
  sig <- res %>%
    filter(padj < padj_cut) %>%
    mutate(
      padj_txt = ifelse(padj < 0.001,
                        formatC(padj, format = "e", digits = 1),   # e.g. 2.3e-05
                        formatC(padj, format = "f", digits = 3)),  # e.g. 0.012
      label = paste0(gsub("_", " ", sub("HALLMARK_", "", pathway)),
                     "  (FDR = ", padj_txt, ")")
    )
  
  if (nrow(sig) == 0) {
    message("Factor ", k, ": no pathways with padj < ", padj_cut, ", skipped")
    return(NULL)
  }
  
  # gene weights per significant pathway
  path_df <- bind_rows(lapply(seq_len(nrow(sig)), function(i) {
    genes <- intersect(pathways[[sig$pathway[i]]], names(stats))
    tibble(label  = sig$label[i],
           NES    = sig$NES[i],
           weight = unname(stats[genes]))
  }))
  
  # background: all genes (NES = NA -> grey)
  bg_df <- tibble(label = "All genes", NES = NA_real_, weight = unname(stats))
  
  # order: highest NES at the top, "All genes" at the bottom
  lvls <- c("All genes", sig$label[order(sig$NES)])
  plot_df <- bind_rows(path_df, bg_df) %>%
    mutate(label = factor(label, levels = lvls))
  
  ggplot(plot_df, aes(x = weight, y = label, fill = NES)) +
    geom_boxplot(orientation = "y", alpha = 0.85, colour = "grey30",
                 linewidth = 0.4, width = 0.7,
                 outlier.size = 0.6, outlier.alpha = 0.5) +
    scale_fill_gradientn(colours = c("#FFF3B0", "#FDB863", "#E66101", "#B2182B"),
                         na.value = "grey80", name = "NES",
                         limits = nes_limits) +
    labs(title = paste0("Factor ", k, " - weight distribution of enriched pathways"),
         x = "Absolute weight", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 9))
}

box_list <- lapply(box_factors, make_box_plot)
names(box_list) <- paste0("Factor", box_factors)
box_list <- Filter(Negate(is.null), box_list)

# view
for (nm in names(box_list)) print(box_list[[nm]])

# save each to its own file (height scales with the number of pathways)
for (nm in names(box_list)) {
  n_rows <- length(unique(box_list[[nm]]$data$label))
  ggsave(paste0("boxplot_fgsea_", nm, ".png"), box_list[[nm]],
         width = 8, height = max(4, 0.4 * n_rows + 1.5), dpi = 300, bg = "white")
}


################################################################################
#check what happens to pathways if directionality is added
# ===== Signed (non-absolute) FGSEA =====
# 1. Signed ranked stats (rna_weight itself was never overwritten) -----------
signed_df <- as.data.frame(rna_weight)
signed_df$gene <- sub("_view_1", "", rownames(signed_df))

signed_factors <- c(1, 8, 9, 11, 12, 13, 14, 15)
signed_stats <- lapply(signed_factors, function(k)
  setNames(signed_df[[paste0("Factor", k)]], signed_df$gene))
names(signed_stats) <- paste0("Factor", signed_factors)

# 2. FGSEA on signed weights -------------------------------------------------
# "std" = two-sided: positive NES = set enriched at the positive end,
#                    negative NES = set enriched at the negative end
set.seed(42)
signed_fgsea <- lapply(signed_stats, function(s)
  fgsea(pathways = pathways, stats = s,
        minSize = 10, maxSize = 500, scoreType = "std"))

# Optional: tails separately (can reveal sets hidden by cancellation in "std")
signed_fgsea_pos <- lapply(signed_stats, function(s)
  fgsea(pathways = pathways, stats = s, minSize = 10, maxSize = 500, scoreType = "pos"))
signed_fgsea_neg <- lapply(signed_stats, function(s)
  fgsea(pathways = pathways, stats = s, minSize = 10, maxSize = 500, scoreType = "neg"))

# Largest |NES| among significant pathways across all factors -> shared symmetric limit
global_nes_lim <- ceiling(max(unlist(lapply(signed_fgsea, function(r)
  abs(r$NES[!is.na(r$padj) & r$padj < 0.05])))))
global_nes_lim   # check the value

# 3. Significant pathways for Factor 11 (sign of NES = direction) -------------
signed_fgsea$Factor11 %>%
  filter(padj < 0.05) %>%
  arrange(desc(NES)) %>%
  select(pathway, NES, padj, size)

# Same for the tail-specific runs
signed_fgsea_pos$Factor11 %>% filter(padj < 0.05) %>% arrange(padj) %>% select(pathway, NES, padj)
signed_fgsea_neg$Factor11 %>% filter(padj < 0.05) %>% arrange(padj) %>% select(pathway, NES, padj)

# 4. Are the TGF-beta genes coherent in sign? --------------------------------
tgfb_genes <- intersect(pathways[["HALLMARK_TGF_BETA_SIGNALING"]], names(signed_stats$Factor11))
tgfb_w <- signed_stats$Factor11[tgfb_genes]
table(sign = ifelse(tgfb_w > 0, "positive", "negative"))   # split of all set members

# leading-edge genes from the std run, with their signed weights
tgfb_le <- signed_fgsea$Factor11 %>%
  filter(pathway == "HALLMARK_TGF_BETA_SIGNALING") %>%
  pull(leadingEdge) %>% unlist()
signed_tgfb_le_df <- tibble(gene = tgfb_le, weight = signed_stats$Factor11[tgfb_le]) %>%
  arrange(desc(weight))
print(signed_tgfb_le_df, n = Inf)

# Same check for the interferon sets, to see whether they share a sign with TGF-beta
for (p in c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE")) {
  w <- signed_stats$Factor11[intersect(pathways[[p]], names(signed_stats$Factor11))]
  cat(p, ": ", sum(w > 0), " positive, ", sum(w < 0), " negative\n", sep = "")
}

# 5. Signed ridge plot --------------------------------------------------------
make_signed_ridge <- function(k, padj_cut = 0.05) {
  res   <- signed_fgsea[[paste0("Factor", k)]]
  stats <- signed_stats[[paste0("Factor", k)]]
  
  sig <- res %>%
    filter(padj < padj_cut) %>%
    mutate(padj_txt = ifelse(padj < 0.001,
                             formatC(padj, format = "e", digits = 1),
                             formatC(padj, format = "f", digits = 3)),
           label = paste0(gsub("_", " ", sub("HALLMARK_", "", pathway)),
                          "  (FDR = ", padj_txt, ")"))
  if (nrow(sig) == 0) { message("Factor ", k, ": nothing significant"); return(NULL) }
  
  path_df <- bind_rows(lapply(seq_len(nrow(sig)), function(i) {
    genes <- intersect(pathways[[sig$pathway[i]]], names(stats))
    tibble(label = sig$label[i], NES = sig$NES[i], weight = unname(stats[genes]))
  }))
  bg_df <- tibble(label = "All genes", NES = NA_real_, weight = unname(stats))
  
  lvls <- c("All genes", sig$label[order(sig$NES)])
  plot_df <- bind_rows(path_df, bg_df) %>% mutate(label = factor(label, levels = lvls))
  
  ggplot(plot_df, aes(x = weight, y = label, fill = NES)) +
    geom_density_ridges(scale = 1.3, alpha = 0.85, colour = "grey30",
                        linewidth = 0.3, trim = TRUE,
                        quantile_lines = TRUE, quantiles = 2) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_fill_gradientn(
      colours = c("#2166AC", "#4393C3", "#92C5DE",   # negative NES: blues
                  "#FFF3B0",                         # NES = 0: pale yellow
                  "#FDB863", "#E66101", "#B2182B"),  # positive NES: your original ramp
      values  = seq(0, 1, length.out = 7),           # evenly spaced; middle colour sits at 0
      limits  = c(-global_nes_lim, global_nes_lim),  # symmetric, same for every plot
      na.value = "grey80", name = "NES") +
    labs(title = paste0("Factor ", k, " - signed weight distribution of enriched pathways"),
         x = "Weight (signed)", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 9))
}

#signed_ridge_list <- Filter(Negate(is.null), lapply(signed_factors, make_signed_ridge))
#names(signed_ridge_list) <- paste0("Factor", signed_factors)[
#  signed_factors %in% signed_factors][seq_along(signed_ridge_list)]


#print(make_signed_ridge(1))
#print(make_signed_ridge(8))
#print(make_signed_ridge(9))
#print(make_signed_ridge(11))
#print(make_signed_ridge(12))
#print(make_signed_ridge(13))
#print(make_signed_ridge(14))
#print(make_signed_ridge(15))
######################################################################################


#merge all of these (later integrate proteomic & metabolomic factors)
fg1<- fgsea_results1 %>% filter(padj <0.05) %>% select(pathway, padj, size)
fg1$Factor <- 'Factor 1'
fg7<- fgsea_results7 %>% filter(padj <0.05) %>% select(pathway, padj, size)
fg7$Factor <- 'Factor 7'
fg8<- fgsea_results8 %>% filter(padj <0.05) %>% select(pathway, padj, size)
fg8$Factor <- 'Factor 8'
fg9<- fgsea_results9 %>% filter(padj <0.05) %>% select(pathway, padj, size)
fg9$Factor <- 'Factor 9'
fg11<- fgsea_results11 %>% filter(padj <0.05) %>% select(pathway, padj, size)
fg11$Factor <- 'Factor 11'
fg12<- fgsea_results12 %>% filter(padj <0.05) %>% select(pathway, padj, size)
fg12$Factor <- 'Factor 12'
fg13<- fgsea_results13 %>% filter(padj <0.05) %>% select(pathway, padj, size)
fg13$Factor <- 'Factor 13'
fg14<- fgsea_results14 %>% filter(padj <0.05) %>% select(pathway, padj, size)
fg14$Factor <- 'Factor 14'
fg15<- fgsea_results15 %>% filter(padj <0.05) %>% select(pathway, padj, size)
fg15$Factor <- 'Factor 15'

pathways <- rbind(fg1, fg8,fg9,fg11,fg12,fg13,fg14,fg15)
pathways$modality <- 'RNA-seq'

pathways$pathway <- as.character(pathways$pathway)

#name pathways more easily
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_ADIPOGENESIS', 'Adipogenesis', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_CHOLESTEROL_HOMEOSTASIS', 'Cholesterol homeostasis', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_E2F_TARGETS', 'E2F targets', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_GLYCOLYSIS', 'Glycolysis', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_INFLAMMATORY_RESPONSE', 'Inflammatory response', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_KRAS_SIGNALING_UP', 'Increased KRAS signalling', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_OXIDATIVE_PHOSPHORYLATION', 'Oxidative phosphorylation', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_TGF_BETA_SIGNALING', 'TGF-B signalling', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_ALLOGRAFT_REJECTION', 'Allograft rejection', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_COAGULATION', 'Coagulation', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', 'Epithelial Mesenchymal Transition', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_HEME_METABOLISM', 'Haem metabolism', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_INTERFERON_ALPHA_RESPONSE', 'IFNa response', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_MITOTIC_SPINDLE', 'Mitotic spindle', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_P53_PATHWAY', 'P53 pathway', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_TNFA_SIGNALING_VIA_NFKB', 'TNFa signalling via NFKB', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_APICAL_JUNCTION', 'Apical junction', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_COMPLEMENT', 'Complement', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_ESTROGEN_RESPONSE_LATE', 'Late estrogen response', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_HYPOXIA', 'Hypoxia', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_INTERFERON_GAMMA_RESPONSE', 'IFNy response', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_MTORC1_SIGNALING', 'MTORC1 signalling', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_PI3K_AKT_MTOR_SIGNALING', 'PI3, AKT, MTOR signalling', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_UV_RESPONSE_DN', 'UV response DN', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_BILE_ACID_METABOLISM', 'Bile acid metabolism', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_DNA_REPAIR', 'DNA repair', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_G2M_CHECKPOINT', 'G2M checkpoint', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_IL6_JAK_STAT3_SIGNALING', 'IL6-JAK-STAT3 signalling', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_KRAS_SIGNALING_DN', 'KRAS signalling DN', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_MYOGENESIS', 'Myogenesis', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_SPERMATOGENESIS', 'Spermatogenesis', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_UV_RESPONSE_UP', 'Increased UV response', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_ANDROGEN_RESPONSE', 'Androgen response', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_APICAL_SURFACE', 'Apical surface', pathways$pathway)
pathways$pathway <- ifelse(pathways$pathway == 'HALLMARK_ESTROGEN_RESPONSE_EARLY', 'Early estrogen response', pathways$pathway)

#note think about a loop for this later on! Must be more efficient, but also need to
#look at the details

#======================================================
#manually adding pathway analysis from below to combine
#note: the actual pathway extraction for the pathways happens on STRING! 

#NOTE HERE: for AAB, given random (lowest) padj value, as there is no overrepresentation analysis, but just taken the most defining AAB
#this needs to be noted in the figure legend

# Create dataframe
df_forbind <- data.frame(
  Factor = c("Factor 2", "Factor 3", "Factor 4", "Factor 5", "Factor 6", "Factor 10", "Factor 7"),
  
  pathway = c("Response to stimulus", "Wnt signalling", "Metabolic process", "Positive regulation of immune response", "Response to endogenous stimulus",  "MYC targets V1", "Histone AABs"),
  
  padj = as.numeric(gsub(",", ".", c("3,0E-08", "0,0174", "4,3E-10", "1,0E-05", "6,8E-03", "0,00019", "0.04"))),
  
  modality = c("Proteomics", "Proteomics", "Proteomics", "Proteomics", "Proteomics", "Proteomics", "AAB+"),
  
  size = NA  # Since "size" column is empty in the image
)

#=========================================================
# RJ - this didnt run

pathways <- rbind(pathways, df_forbind)

#now select 5 most significant pathways per factor (if they pass threshold of course)

#pathways <- pathways %>% group_by(Factor) %>% arrange(padj, .by_group = TRUE) %>% slice_head(n = 5) %>% ungroup()

top_pathways <- pathways %>%
  group_by(Factor) %>%
  arrange(padj, .by_group = TRUE) %>%
  filter(
    row_number() <= if_else(Factor[1] %in% c('Factor 11', 'Factor 13'), 2L, 1L)
  ) %>%
  ungroup()

# Combine and remove duplicates
pathways <- top_pathways



# Compute -log10(padj)
library(ggrepel)
pathways$log10padj <- -log10(pathways$padj)

#Reorder Y-axis
pathways$Factor <- factor(pathways$Factor, levels = rev(unique(pathways$Factor)))

max_per_factor <- pathways %>%
  group_by(Factor) %>%
  summarise(max_x = max(log10padj), .groups = 'drop')

#Merge back for plotting
pathways <- left_join(pathways, max_per_factor, by = "Factor")

#Threshold line at padj = 0.05
sig_threshold <- -log10(0.05)

#relevel factors
ordered_levels <- paste("Factor", 15:1)
pathways$Factor <- factor(pathways$Factor, levels = ordered_levels)


# Plot
ppath <- ggplot(pathways, aes(x = log10padj, y = Factor)) +
  # Horizontal lines per factor up to max x (no need to convert y to numeric)
  geom_segment(
    aes(x = 0, xend = max_x, y = Factor, yend = Factor),
    color = "black", alpha = 0.2, linewidth = 0.5
  ) +
  
  # Dots with fixed size and colored by modality
  geom_point(aes(color = modality), size = 5, alpha = 0.8) +
  
  # Pathway labels
  geom_text_repel(aes(label = pathway),
                  size = 4,
                  max.overlaps = Inf,     
                  box.padding = 0.5,      
                  point.padding = 0.5,    
                  force = 3,              
                  force_pull = 1,         
                  segment.color = "grey50",
                  segment.size = 0.3,
                  min.segment.length = 0) +
  # Vertical line for significance
  geom_vline(xintercept = sig_threshold, linetype = "dotted", color = "red", linewidth = 0.8) +
  
  labs(
    x = expression(-log[10](FDR~adjusted~p~value)),
    y = "Factor",
    title = "Enriched Pathways by Factor",
    color = "Modality"
  ) +
  theme_pubr() +
  theme(
    axis.text.y = element_text(size = 13),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = c(0.8, 0.5), legend.justification = c("left", "center"))


ggarrange(ppath, survp2, ncol = 2, widths = c(1.3, 1), align = "hv") +
  theme(plot.margin = margin(0, 0, 0, 0))

ppath <- ppath + theme(plot.margin = margin(t = 5, r = 2, b = 5, l = 5))
survp2 <- survp2 + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 2))
ggarrange(ppath, survp2, ncol = 2, widths = c(1.3, 1), align = "hv")

library(patchwork)
ppath <- ppath + theme(plot.margin = margin(t = 5, r = 0, b = 5, l = 5))
survp2 <- survp2 + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 0))

combined <- ppath + survp2 + 
  plot_layout(ncol = 2, widths = c(1.3, 1))
combined
#=================================================================================================================
#EdB continued on 23/10/24
#investigate associations between factors and drugs/combidities
clean_comorbidity <- readRDS("/mnt/session_data/mounted-data-readonly/PAH-OpenClinica/20240726-freeze/clean_comorbidity.rds")
clean_medication <- readRDS("/mnt/session_data/mounted-data-readonly/PAH-OpenClinica/20240726-freeze/clean_medication.rds")

#do analysis of factor association with for comorbid disorders
comorbid_wide <- clean_comorbidity %>% distinct(id, cfh_comorbid_disease_diagnosis) %>% mutate(value = "yes") %>%  pivot_wider(names_from = cfh_comorbid_disease_diagnosis, values_from = value, values_fill = "no") 

#some patients do not seem to have any data, these do not have any comorbidities 
def_scores <- definitive_mofa_model_factor_scores %>% filter(!samples %in% comorbid_wide$id)
#these 18 patients need to be added
add <- def_scores %>% select(samples)
add <- add %>% rename(id = samples)
comorbid_wide <- bind_rows(comorbid_wide, add) %>%
  mutate(across(-id, ~ replace_na(.x, "no")))

c1 <- merge(comorbid_wide, factor_scores_definitive, by.x='id', by.y='samples')

#if entire col is no for med, remove it (no use in testing it)
c1 <- c1 %>% select_if(~ !all(. == "no"))

#also needs to have >5 yes' because otherwise not informative
# RJ - do you not think 5 is a bit low?
columns_to_keep <- sapply(2:293, function(i) {
  sum(c1[[i]] == "yes") > 8
})

# Select those columns along with the numeric columns (302:316)
c1 <- c1 %>% select(1, which(columns_to_keep) + 1, 294:308)

age_sex <- clin %>% select(id, age_diagnosis, sex)
c1 <- merge(c1, age_sex, by='id')


med_cols <- colnames(c1)[2:17]
factor_cols <- colnames(c1)[18:32]

# Function to run an adjusted linear model and extract the p-value
run_adjusted_model <- function(char_col, numeric_cols, data) {
  sapply(numeric_cols, function(num_col) {
    tryCatch({
      # Construct formula dynamically: Factor ~ Medication + age + sex
      # Backticks (` `) are used to safely handle any spaces or special characters in column names
      form <- as.formula(paste0("`", num_col, "` ~ `", char_col, "` + age_diagnosis + sex"))
      
      # Fit the linear model
      fit <- lm(form, data = data)
      
      # Extract the coefficient matrix
      coefs <- summary(fit)$coefficients
      
      # Find the row corresponding to the medication's effect (e.g., 'loratidineyes')
      # We use grep to find the term that contains the medication's column name
      term_idx <- grep(char_col, rownames(coefs))
      
      if (length(term_idx) > 0) {
        # Return the p-value (which is in the 4th column of the summary matrix)
        return(coefs[term_idx[1], 4])
      } else {
        return(NA)
      }
    }, error = function(e) {
      return(NA) # Return NA if the model fails (e.g., missing data or singular fit)
    })
  })
}

# Apply the function across all medication columns
pvals_adjusted <- lapply(med_cols, function(med) {
  run_adjusted_model(med, factor_cols, data = c1)
})

# Convert the resulting list into a dataframe
pval_df <- do.call(rbind, pvals_adjusted)

# Set the row names (medications) and column names (factors)
rownames(pval_df) <- med_cols
colnames(pval_df) <- factor_cols

pval_df <- as.data.frame(pval_df)

# Apply FDR correction to the columns
pval_df_fdr <- as.data.frame(apply(pval_df, 2, function(col) {
  p.adjust(col, method = "fdr")
}))

#########################################################################
#EdB add: heatmap with significance:

# Number of patients with each comorbidity (count of "yes" in c1)
n_yes <- sapply(med_cols, function(m) sum(c1[[m]] == "yes", na.rm = TRUE))
n_df  <- data.frame(Comorbidity_raw = med_cols, n = n_yes)

# Helper: keep only the number after "Factor" (e.g. "filtered__group1.Factor12" -> "12")
factor_num <- function(x) sub(".*Factor", "", x)

#make names nicer
pretty_map <- c(
  "gerd" = "GERD", "dm2" = "DM2", "osa" = "OSA", "ckd" = "CKD",
  "copd" = "COPD", "pfo" = "PFO", "afib" = "AFib",
  "vit d deficiency" = "Vitamin D deficiency"
)

nice_name <- function(x) {
  x   <- tolower(gsub("_", " ", x))
  out <- ifelse(x %in% names(pretty_map), pretty_map[x],
                paste0(toupper(substring(x, 1, 1)), substring(x, 2)))
  unname(out)
}

# Long format
fdr_long <- pval_df_fdr %>%
  mutate(Comorbidity_raw = rownames(pval_df_fdr)) %>%
  pivot_longer(-Comorbidity_raw, names_to = "Factor_raw", values_to = "FDR") %>%
  left_join(n_df, by = "Comorbidity_raw") %>%
  mutate(
    neglog = -log10(FDR),
    stars = case_when(
      is.na(FDR)  ~ "",
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE        ~ ""
    ),
    Comorbidity = paste0(nice_name(Comorbidity_raw), " (n = ", n, ")"),
    Factor = paste("Factor", factor_num(Factor_raw))
  )

# Order comorbidities by their strongest association
row_order <- fdr_long %>%
  group_by(Comorbidity) %>%
  summarise(min_fdr = min(FDR, na.rm = TRUE)) %>%
  arrange(desc(min_fdr)) %>%
  pull(Comorbidity)

# Order factors numerically: Factor 1 ... Factor 15
factor_levels <- paste("Factor", factor_num(factor_cols))

fdr_long <- fdr_long %>%
  mutate(
    Comorbidity = factor(Comorbidity, levels = row_order),
    Factor = factor(Factor, levels = factor_levels)
  )

heat_fdr <- ggplot(fdr_long, aes(x = Factor, y = Comorbidity, fill = neglog)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = stars), size = 5, vjust = 0.75, colour = "black") +
  scale_fill_gradient(low = "white", high = "#B2182B",
                      na.value = "grey90",
                      name = expression(-log[10](FDR))) +
  labs(x = NULL, y = NULL,
       title = "Association of factors with comorbidities",
       subtitle = "Adjusted for age and sex. * FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

print(heat_fdr)

ggsave("heatmap_factor_comorbidity_FDR.png", heat_fdr,
       width = 9, height = max(4, 0.3 * nlevels(fdr_long$Comorbidity) + 2),
       dpi = 300, bg = "white")

###########################################################################
# Convert the result back to a dataframe
pval_df_fdr <- as.data.frame(pval_df_fdr)
write.csv2(pval_df_fdr, file='FDR_adj_pvals_lm_comorbs_for_factors_9March26.csv')

#NOTE HERE: factor 15 seems to have become a thyroid factor (?!)

#now plot
p1 <- ggplot(c1, aes(x = DM2, y = filtered__group1.Factor2, fill = DM2)) +
  geom_boxplot() +
  scale_fill_manual(values = c("yes" = "blue", "no" = "red")) +
  labs(x = "T2DM", y = "Factor 2") +
  #stat_compare_means(method='t.test', label.y=3) +
  ggtitle('Factor 2 for T2DM, FDR = 0.023') +
  theme_pubr()

p2 <- ggplot(c1, aes(x = GERD, y = filtered__group1.Factor4, fill = GERD)) +
  geom_boxplot() +
  scale_fill_manual(values = c("yes" = "blue", "no" = "red")) +
  labs(x = "GERD", y = "Factor 4") +
  #stat_compare_means(method='t.test', label.y=3) +
  ggtitle('Factor 4 for GERD, FDR <0.0001') +
  theme_pubr()

p3 <- ggplot(c1, aes(x = obesity, y = filtered__group1.Factor7, fill = obesity)) +
  geom_boxplot() +
  scale_fill_manual(values = c("yes" = "blue", "no" = "red")) +
  labs(x = "Obesity", y = "Factor 7") +
  #stat_compare_means(method='t.test', label.y=3) +
  ggtitle('Factor 7 for obesity, FDR = 0.049') +
  theme_pubr()


p4 <- ggplot(c1, aes(x = hypothyroidism, y = filtered__group1.Factor15, fill = hypothyroidism)) +
  geom_boxplot() +
  scale_fill_manual(values = c("yes" = "blue", "no" = "red")) +
  labs(x = "Hypothyroidism", y = "Factor 15") +
  #stat_compare_means(method='t.test', label.y=1.5) +
  ggtitle('Factor 15 for hypothyroidism, FDR = 0.012') +
  theme_pubr()

g1 <- ggarrange(p1, p2,p3,p4, nrow = 2, ncol=2)


#RJ - i think all this data would look great together in a heatmap - very rough version!
log_pval_df <- -log10(pval_df_fdr)
mat <- as.matrix(log_pval_df)

pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  border_color = NA
)

#=========================================================
#repeat this for drug classes
clean_medication_wide <- clean_medication %>% distinct(id, class) %>% mutate(value = "yes") %>%  pivot_wider(names_from = class, values_from = value, values_fill = "no") 
c2 <- merge(clean_medication_wide, factor_scores_definitive, by.x='id', by.y='samples')

#remove the NA col
c2 <- c2[,-3]

#if entire col is no for med, remove it (no use in testing it)
c2 <- c2 %>% select_if(~ !all(. == "no"))


#also needs to have >8 yes' because otherwise not informative
columns_to_keep <- sapply(2:49, function(i) {
  sum(c2[[i]] == "yes") > 8
})

# Select those columns along with the numeric columns (302:316)
c2 <- c2 %>% select(1, which(columns_to_keep) + 1, 50:64)
age_sex <- clin %>% select(id, age_diagnosis, sex)
c2 <- merge(c2, age_sex, by='id')

#note, most look roughly normal! 

#ADJUST FOR AGE & SEX
# Define your columns based on your dataset
# (Double-check these indices after adding age and sex!)
med_cols <- colnames(c2)[2:33]
factor_cols <- colnames(c2)[34:48]

# RJ - are you sure this has adjusted for age and sex? I find these functions quite confusing and this bit - "` + age_diagnosis + sex" makes me worried its not right, did you compare without adjustment to show its worked?
#EdB comment: double checked and this is correct 

# Function to run an adjusted linear model and extract the p-value
run_adjusted_model <- function(char_col, numeric_cols, data) {
  sapply(numeric_cols, function(num_col) {
    tryCatch({
      # Construct formula dynamically: Factor ~ Medication + age + sex
      # Backticks (` `) are used to safely handle any spaces or special characters in column names
      form <- as.formula(paste0("`", num_col, "` ~ `", char_col, "` + age_diagnosis + sex"))
      
      # Fit the linear model
      fit <- lm(form, data = data)
      
      # Extract the coefficient matrix
      coefs <- summary(fit)$coefficients
      
      # Find the row corresponding to the medication's effect (e.g., 'loratidineyes')
      # We use grep to find the term that contains the medication's column name
      term_idx <- grep(char_col, rownames(coefs))
      
      if (length(term_idx) > 0) {
        # Return the p-value (which is in the 4th column of the summary matrix)
        return(coefs[term_idx[1], 4])
      } else {
        return(NA)
      }
    }, error = function(e) {
      return(NA) # Return NA if the model fails (e.g., missing data or singular fit)
    })
  })
}

# Apply the function across all medication columns
pvals_adjusted <- lapply(med_cols, function(med) {
  run_adjusted_model(med, factor_cols, data = c2)
})

# Convert the resulting list into a dataframe
pval_df <- do.call(rbind, pvals_adjusted)

# Set the row names (medications) and column names (factors)
rownames(pval_df) <- med_cols
colnames(pval_df) <- factor_cols

pval_df <- as.data.frame(pval_df)

# Apply FDR correction to the columns
pval_df_fdr <- as.data.frame(apply(pval_df, 2, function(col) {
  p.adjust(col, method = "fdr")
}))

####################################################################
#EdB add medication heatmap
# Number of patients on each drug class (count of "yes" in c2)
n_yes_med <- sapply(med_cols, function(m) sum(c2[[m]] == "yes", na.rm = TRUE))
n_df_med  <- data.frame(Medication_raw = med_cols, n = n_yes_med)

# Helper: keep only the number after "Factor"
factor_num <- function(x) sub(".*Factor", "", x)

# Make medication names nicer (keys are lowercase versions of the raw names)
pretty_map_med <- c(
  "su derrivate"             = "Sulfonylurea",
  "ca channel blocker"       = "Calcium channel blocker",
  "h2 blocker"               = "H2 blocker",
  "glp1 agonist"             = "GLP-1 agonist",
  "ssri snri"                = "SSRI/SNRI",
  "acei arbs"                = "ACEi/ARBs",
  "longterm oxygen"          = "Long-term oxygen",
  "anti platelet"            = "Antiplatelet",
  "heparin lmwh"             = "Heparin (LMWH)",
  "immunosuppresant other"   = "Other immunosuppressant",
  "iron suppletion"          = "Iron supplementation",
  "pain other"               = "Other analgesic",
  "imantinib"                = "Imatinib",
  "benzo"                    = "Benzodiazepine",
  "pde5i"                    = "PDE5i",
  "era"                      = "ERA",
  "pra"                      = "PRA",
  "sgc"                      = "sGC",
  "tca"                      = "TCA",
  "ppi"                      = "PPI",
  "nsaid"                    = "NSAID"
)

nice_name_med <- function(x) {
  x   <- gsub("_", " ", x)
  key <- tolower(x)
  out <- ifelse(key %in% names(pretty_map_med), pretty_map_med[key],
                paste0(toupper(substring(x, 1, 1)), substring(x, 2)))
  unname(out)
}

# Long format
fdr_long_med <- pval_df_fdr %>%
  mutate(Medication_raw = rownames(pval_df_fdr)) %>%
  pivot_longer(-Medication_raw, names_to = "Factor_raw", values_to = "FDR") %>%
  left_join(n_df_med, by = "Medication_raw") %>%
  mutate(
    neglog = -log10(FDR),
    stars = case_when(
      is.na(FDR)  ~ "",
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE        ~ ""
    ),
    Medication = paste0(nice_name_med(Medication_raw), " (n = ", n, ")"),
    Factor = paste("Factor", factor_num(Factor_raw))
  )

# Order medications by their strongest association
row_order_med <- fdr_long_med %>%
  group_by(Medication) %>%
  summarise(min_fdr = min(FDR, na.rm = TRUE)) %>%
  arrange(desc(min_fdr)) %>%
  pull(Medication)

# Order factors numerically: Factor 1 ... Factor 15
factor_levels <- paste("Factor", factor_num(factor_cols))

fdr_long_med <- fdr_long_med %>%
  mutate(
    Medication = factor(Medication, levels = row_order_med),
    Factor = factor(Factor, levels = factor_levels)
  )

heat_fdr_med <- ggplot(fdr_long_med, aes(x = Factor, y = Medication, fill = neglog)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = stars), size = 5, vjust = 0.75, colour = "black") +
  scale_fill_gradient(low = "white", high = "#B2182B",
                      na.value = "grey90",
                      name = expression(-log[10](FDR))) +
  labs(x = NULL, y = NULL,
       title = "Association of factors with medication classes",
       subtitle = "Adjusted for age and sex. * FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

print(heat_fdr_med)

ggsave("heatmap_factor_medication_FDR.png", heat_fdr_med,
       width = 9, height = max(4, 0.3 * nlevels(fdr_long_med$Medication) + 2),
       dpi = 300, bg = "white")

ggarrange(heat_fdr, heat_fdr_med)
#####################################################################

write.csv2(pval_df_fdr, file='FDR_adj_pvals_lm_med_classes_for_factors_9March26.csv')


p1 <- ggplot(c2, aes(x = Diuretics, y = filtered__group1.Factor3, fill = Diuretics)) +
  geom_boxplot() +
  scale_fill_manual(values = c("yes" = "blue", "no" = "red")) +
  labs(x = "Diuretics", y = "Factor 3") +
  ggtitle('Factor 3 for diuretics usage, FDR = 0.000010') +
  theme_pubr()

p2 <- ggplot(c2, aes(x = PA, y = filtered__group1.Factor3, fill = PA)) +
  geom_boxplot() +
  scale_fill_manual(values = c("yes" = "blue", "no" = "red")) +
  labs(x = "Prostacyclin annalogues", y = "Factor 3") +
  ggtitle('Factor 3 for Prostacyclin annalogue medication, FDR = 0.025') +
  theme_pubr()

p3 <- ggplot(c2, aes(x = topical_corticosteroid, y = filtered__group1.Factor4, fill = topical_corticosteroid)) +
  geom_boxplot() +
  scale_fill_manual(values = c("yes" = "blue", "no" = "red")) +
  labs(x = "Topical steroids", y = "Factor 4") +
  ggtitle('Factor 4 for topical steroid medication, FDR = 0.0021') +
  theme_pubr()

p4 <- ggplot(c2, aes(x = antacid, y = filtered__group1.Factor5, fill = antacid)) +
  geom_boxplot() +
  scale_fill_manual(values = c("yes" = "blue", "no" = "red")) +
  labs(x = "Hypothyroid medication usage", y = "Factor 5") +
  ggtitle('Factor 5 for antacid medication, FDR = 0.010') +
  theme_pubr()

p5 <- ggplot(c2, aes(x = SU_derrivate, y = filtered__group1.Factor7, fill = SU_derrivate)) +
  geom_boxplot() +
  scale_fill_manual(values = c("yes" = "blue", "no" = "red")) +
  labs(x = "Sulfonylurea medication usage", y = "Factor 7") +
  ggtitle('Factor 7 for Sulfonylurea medication, FDR = 0.0014') +
  theme_pubr()

p6 <- ggplot(c2, aes(x = hypothyroid_drugs, y = filtered__group1.Factor15, fill = hypothyroid_drugs)) +
  geom_boxplot() +
  scale_fill_manual(values = c("yes" = "blue", "no" = "red")) +
  labs(x = "Hypothyroid drug usage", y = "Factor 15") +
  ggtitle('Factor 15 for hypothyroidism medication, FDR = 0.026') +
  theme_pubr()


ggarrange(p1,p2,p3,p4,p5,p6, ncol=3, nrow=2)

log_pval_df <- -log10(pval_df_fdr)
mat <- as.matrix(log_pval_df)

pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  border_color = NA
)

#==================================================================================================
#STRING analysis
#select proteins for selection
pw <- proteome_weight %>% select(Factor2, Factor3, Factor4, Factor5, Factor6, Factor7, Factor10)
pw$protein <- rownames(pw)
getwd()

#remove _view_2 from this
pw[] <- lapply(pw, function(x) if(is.character(x)) sub("_view_2", "", x) else x)
proteomics_proteins_to_keep <- readRDS("/mnt/session_data/mounted-data-readonly/proteomics_proteins_to_keep.rds")
#saveRDS(pw, 'proteomics_factor_scores_definitive_MOFA_for_String.rds')

pw2 <- merge(pw, proteomics_proteins_to_keep, by.x='protein', by.y='id_protein')

#now select those within filter
#this filter is arbitrary
# RJ - looking at the histograms and the number of proteins here, I am concerned this is too stringent. Can we not change it to the top 150 most impactful proteins (top 75/bottom 75) ?

pw2_factor2 <- pw2 %>% filter(Factor2 >0.45 | Factor2 < -0.45)
pw2_factor2 <- pw2_factor2 %>% select(protein, Factor2, EntrezGeneSymbol)
pw2_factor2$Factor2 <- 'yes'

pw2_factor3 <- pw2 %>% filter(Factor3 >0.45 | Factor3 < -0.45)
pw2_factor3 <- pw2_factor3 %>% select(protein, Factor3, EntrezGeneSymbol)
pw2_factor3$Factor3 <- 'yes'

pw2_factor4 <- pw2 %>% filter(Factor4 >0.45 | Factor4 < -0.45)
pw2_factor4 <- pw2_factor4 %>% select(protein, Factor4, EntrezGeneSymbol)
pw2_factor4$Factor4 <- 'yes'

pw2_factor5 <- pw2 %>% filter(Factor5 >0.45 | Factor5 < -0.45)
pw2_factor5 <- pw2_factor5 %>% select(protein, Factor5, EntrezGeneSymbol)
pw2_factor5$Factor5 <- 'yes'

pw2_factor6 <- pw2 %>% filter(Factor6 >0.45 | Factor6 < -0.45)
pw2_factor6 <- pw2_factor6 %>% select(protein, Factor6, EntrezGeneSymbol)
pw2_factor6$Factor6 <- 'yes'

pw2_factor7 <- pw2 %>% filter(Factor7 >0.45 | Factor7 < -0.45)
pw2_factor7 <- pw2_factor7 %>% select(protein, Factor7, EntrezGeneSymbol)
#pw2_factor7$Factor7 <- 'yes'
#noting retained


pw2_factor10 <- pw2 %>% filter(Factor10 >0.45 | Factor10 < -0.45)
pw2_factor10 <- pw2_factor10 %>% select(protein, Factor10, EntrezGeneSymbol)
pw2_factor10$Factor10 <- 'yes'
#NOTE: for this taking the top 3 Biological processes or KEGG pathways

# not sure why youve merged these
prots <- merge(pw2_factor2, pw2_factor3, by=c('protein', 'EntrezGeneSymbol'), all=T)
prots2 <- merge(pw2_factor4, pw2_factor5, by=c('protein', 'EntrezGeneSymbol'), all=T)
prots3 <- merge(pw2_factor6, pw2_factor10, by=c('protein', 'EntrezGeneSymbol'), all=T)
protsa <- merge(prots2, prots3, by=c('protein', 'EntrezGeneSymbol'), all=T)
prots <- merge(prots, protsa, by=c('protein', 'EntrezGeneSymbol'), all=T)
write_rds(prots, file='proteins_takn_forward_for_STRING_6Nov25.rds')
write.csv2(prots, file='proteins_takn_forward_for_STRING_6Nov25.csv')
#=========================================================================================
#metabolomic analysis 
metabolites_conversion_file <- readRDS("/mnt/session_data/mounted-data-readonly/metabolites_conversion_file.rds")

metab_weight$id_metabolite <- rownames(metab_weight)

metab_weight <- merge(metab_weight, metabolites_conversion_file, by='id_metabolite')
meta_check <- metab %>% filter(!rownames(metab) %in% metab_weight$id_metabolite)

full_table <- metab_weight %>% count(SUB_PATHWAY)

#these are only unnamed metabolites, so fine for pathway analysis
f2m <- metab_weight %>% filter(abs(Factor2) >0.3)
summary(as.factor(f2m$SUB_PATHWAY))

filtered_table2 <- f2m %>%count(SUB_PATHWAY)

#do overrepresentation analysis
f3m <- metab_weight %>% filter(abs(Factor3) >0.3)
summary(as.factor(f3m$SUB_PATHWAY))

filtered_table3 <- f3m %>% count(SUB_PATHWAY)

f4m <- metab_weight %>% filter(abs(Factor4) >0.3)
summary(as.factor(f4m$SUB_PATHWAY))

filtered_table4 <- f4m %>% count(SUB_PATHWAY)

f5m <- metab_weight %>% filter(abs(Factor5) >0.3)
summary(as.factor(f5m$SUB_PATHWAY))

filtered_table5 <- f5m %>% count(SUB_PATHWAY)

f6m <- metab_weight %>% filter(abs(Factor6) >0.3)
summary(as.factor(f6m$SUB_PATHWAY))

filtered_table6 <- f6m %>% count(SUB_PATHWAY)


f7m <- metab_weight %>% filter(abs(Factor7) >0.3)
summary(as.factor(f7m$SUB_PATHWAY))

filtered_table7 <- f7m %>% count(SUB_PATHWAY)


f10m <- metab_weight %>% filter(abs(Factor10) >0.3)
summary(as.factor(f10m$SUB_PATHWAY))

filtered_table10 <- f10m %>% count(SUB_PATHWAY)


mdf <- metab_weight
mdf$f2 <- ifelse(metab_weight$id_metabolite %in% f2m$id_metabolite, 'YES', 'NO')
mdf$f3 <- ifelse(metab_weight$id_metabolite %in% f3m$id_metabolite, 'YES', 'NO')
mdf$f6 <- ifelse(metab_weight$id_metabolite %in% f6m$id_metabolite, 'YES', 'NO')
mdf$f10 <- ifelse(metab_weight$id_metabolite %in% f10m$id_metabolite, 'YES', 'NO')
mdf$metab_in_pathway <- 'YES'

mdf2 <- mdf %>% pivot_wider(names_from = SUB_PATHWAY, values_from = metab_in_pathway)

#now make NA no
mdf2[is.na(mdf2)] <- "NO"

tt=with(mdf2, table(f2, `Polyamine Metabolism`))
fisher.test(tt)

#this works, now try with loop
fisher_results <- lapply(names(mdf2)[35:84], function(col_name) {
  tt <- with(mdf2, table(f10, get(col_name))) 
  test_result <- fisher.test(tt)  
  
  data.frame(Column = col_name, P_Value = test_result$p.value)
})


fisher_results_df <- do.call(rbind, fisher_results)
#now FDR
fisher_results_df$fdr <- p.adjust(fisher_results_df$P_Value, method='fdr')
write.csv2(fisher_results_df, file='Factor2_metabolite_enrichment_Fisher_test_6Nov25.csv')


# RJ - redoing factor 2

fisher_results <- lapply(names(mdf2)[35:84], function(col_name) {
  tt <- with(mdf2, table(f2, get(col_name))) 
  test_result <- fisher.test(tt)  
  
  data.frame(Column = col_name, P_Value = test_result$p.value)
})


fisher_results_df <- do.call(rbind, fisher_results)
#now FDR
fisher_results_df$fdr <- p.adjust(fisher_results_df$P_Value, method='fdr')
write.csv2(fisher_results_df, file='Factor2_metabolite_enrichment_Fisher_test_6Nov25.csv')

#new loop for factor 3
fisher_results <- lapply(names(mdf2)[35:84], function(col_name) {
  tt <- with(mdf2, table(f3, get(col_name))) 
  test_result <- fisher.test(tt)  
  
  data.frame(Column = col_name, P_Value = test_result$p.value)
})


fisher_results_df <- do.call(rbind, fisher_results)
#now FDR
fisher_results_df$fdr <- p.adjust(fisher_results_df$P_Value, method='fdr')
write.csv2(fisher_results_df, file='Factor3_metabolite_enrichment_Fisher_test_6Nov25.csv')


#now factor 6
fisher_results <- lapply(names(mdf2)[35:84], function(col_name) {
  tt <- with(mdf2, table(f6, get(col_name))) 
  test_result <- fisher.test(tt)  
  
  data.frame(Column = col_name, P_Value = test_result$p.value)
})


fisher_results_df <- do.call(rbind, fisher_results)
#now FDR
fisher_results_df$fdr <- p.adjust(fisher_results_df$P_Value, method='fdr')

write.csv2(fisher_results_df, file='Factor6_metabolite_enrichment_Fisher_test_6Nov25.csv')

#now factor 10
fisher_results <- lapply(names(mdf2)[35:84], function(col_name) {
  tt <- with(mdf2, table(f10, get(col_name))) 
  test_result <- fisher.test(tt)  
  
  data.frame(Column = col_name, P_Value = test_result$p.value)
})


fisher_results_df <- do.call(rbind, fisher_results)
#now FDR
fisher_results_df$fdr <- p.adjust(fisher_results_df$P_Value, method='fdr')

write.csv2(fisher_results_df, file='Factor10_metabolite_enrichment_Fisher_test_6Nov25.csv')

#================================================================================================
#validation analyses
#EdB continued here on 19-12-24
#Defintive MOFA - model 2 - validation
#this is CIBER adjusted including leukocytes with  missingness! 

library(tidyverse)
library(survival)
library(survminer)
library(MOFA2)
library(fgsea)
library(msigdbr)
library(data.table)
library(pheatmap)
library(Publish)
getwd()

# RJ - cleared environment first
########################################################################################
#get mortality data ready
clin_dat_at_diagnosis_clean <- readRDS("/mnt/session_data/mounted-data-readonly/PAH-OpenClinica/20240726-freeze/clin_dat_at_diagnosis_clean.rds")
follow_up_data <- readRDS("/mnt/mounted-data/PAH-OpenClinica/20240726-freeze/follow_up_data.rds")

time_df<- clin_dat_at_diagnosis_clean %>% select(id, diagnosis_date)
fudf <- follow_up_data %>% filter(visit == 1) %>% select(id, visit_date)
time_df <- merge(time_df, fudf, by='id')
time_df$visit_date <- as.Date(time_df$visit_date)
time_df$time_to_sampling <- as.numeric(time_df$visit_date - time_df$diagnosis_date)/365.24
time_df <- time_df[,c(1,3,4)]

mort_df <-  clin_dat_at_diagnosis_clean %>% dplyr::select('id', 'diagnosis_verified', 'sex', 'age_diagnosis', 'bs_bmi', 'diagnosis_date', 'sub_cause', 'sub_date', 'hb_pvr_calc', 'hb_rap_m', 'egfr_mdrd', 'cbt_inflammation_crp_mgpl')
mort_df <- merge(mort_df, time_df, by= 'id')
mort_df$sub_date <- as.Date(mort_df$sub_date)
mort_df$diagnosis_date <- as.Date(mort_df$diagnosis_date)
mort_df <- unique(mort_df)
Cenus_1Dec_2023_UK_PAH_cohort <- readRDS("/mnt/mounted-data/Cenus_1Dec_2023_UK_PAH_cohort.rds")

mort_df <- merge(mort_df, Cenus_1Dec_2023_UK_PAH_cohort, by='id')

mort_df$census_cause <- ifelse(mort_df$census_cause == 'DEATH', 'death', mort_df$census_cause)

#make sure sub_date is complete and format the data
mort_df$sub_date <- as.Date(mort_df$death_or_census_date)
mort_df$sub_date2 <- fifelse(is.na(mort_df$sub_date), ymd("2023-12-01"), mort_df$sub_date)
mort_df$surv_time <- (mort_df$sub_date2 - mort_df$diagnosis_date)
mort_df[c('surv_time', 'day')] <- str_split_fixed(mort_df$surv_time, ' ', 2)
mort_df$surv_time <- as.numeric(mort_df$surv_time)
mort_df$surv_time <- mort_df$surv_time/365.2422
mort_df$event <- ifelse(mort_df$census_cause == 'death', 1, 0)
mort_df$event <- ifelse(is.na(mort_df$event), 0, mort_df$event)
mort_df$diagnosis_verified <- as.factor(mort_df$diagnosis_verified)
mort_df$sex <- as.factor(mort_df$sex)
mort_df <- mort_df %>% filter(!is.na(surv_time))
#remove negative survival times (due to census date)
mort_df <- mort_df %>% filter(surv_time >=0)
mort_df <- unique(mort_df)

#censor at 15 years for cohort
mort_df$event <- ifelse(mort_df$surv_time >15 & mort_df$event == 1, 0, mort_df$event)
mort_df$surv_time <- ifelse(mort_df$surv_time >15, 15, mort_df$surv_time)

mort_df_for_later <- mort_df

clin <- clin_dat_at_diagnosis_clean %>% select(id, age_diagnosis, sex, bs_bmi, diagnosis_verified, functional_class, cbt_card_ntprobnp_ngpl, cbt_inflammation_crp_mgpl, cbt_haem_platelets_x10e9pl, cbt_haem_redcell_dist_pct, cbt_haem_wbc_x10e9pl, cbt_haem_hb_gpl, hb_pvr_calc, hb_rap_m, hb_pap_m, hb_pawp_m, hb_cardiac_output_value_1, hb_cardiac_index_value_1, ep_1_distance_meters, hb_pawp_m, lf_fev1_pc, lf_fvc_pc, lf_kco_pc, cbt_thyr_tsh_mupl)
######################################################################
#load the data
load("/mnt/session_data/MOFA_matrix_complete_including_missigness_after_filtering_30Jan25.RData")
#omit this load below, as this is the old dataframe
#load("/mnt/session_data/MOFA_matrix_complete_including_missigness.RData")
load("/mnt/session_data/MOFA_matrix_list_filtered_pts.RData")
DOUBLE_check_patient_WBC_PC_residuals <- readRDS("/mnt/session_data/DOUBLE_check_patient_WBC_PC_residuals.rds")

a1 <- as.data.frame(filtered_matrix_mofa[[1]])

transcript <- DOUBLE_check_patient_WBC_PC_residuals
transcript1 <- na.omit(t(transcript))
transcript1 <- as.data.frame(transcript1)
proteome <- as.data.frame(mofa_matrix_full[[2]])
proteome1 <- na.omit(t(proteome))
proteome1 <- as.data.frame(proteome1)
metab <- as.data.frame(mofa_matrix_full[[3]])
metab1 <- na.omit(t(metab))
metab1 <- as.data.frame(metab1)
aab <- as.data.frame(mofa_matrix_full[[4]])
aab1 <- na.omit(t(aab))
aab1 <- as.data.frame(aab1)
ciber <- as.data.frame(mofa_matrix_full[[5]])

df1 <- as.data.frame(rownames(transcript1))
df1$transcriptome <- 1
names(df1) <- c('id', 'transcriptome')
df2 <- as.data.frame(rownames(proteome1))
df2$proteome <- 1
names(df2) <- c('id', 'proteome')
df3 <- as.data.frame(rownames(metab1))
df3$metabolome <- 1
names(df3) <- c('id', 'metabolome')
df4 <- as.data.frame(rownames(aab1))
df4$aab <- 1
names(df4) <- c('id', 'aab')

d1 <- merge(df1, df2, by='id', all=T)
d2 <- merge(df3, df4, by='id', all=T)
d3 <- merge(d1, d2, by='id', all=T)

d3$sum_row <- rowSums(d3[,2:5], na.rm=T)
d3 <- d3 %>% filter(sum_row >1)

#467 patients left

rownames(d3) <- d3$id

d4 <- as.data.frame(t(d3))

t1 <- transcript %>% select(all_of(names(d4)))
t1 <- as.data.frame(t(t1))
t1[,1:4524] <- lapply(t1[,1:4524], as.numeric)
t1 <- scale(t1, scale=TRUE, center=TRUE)
t1 <- as.data.frame(t(t1))
t1 <- as.matrix(t1)
p1 <- proteome %>% select(all_of(names(d4)))
p1 <- as.data.frame(t(p1))
p1 <- scale(p1, scale=TRUE, center=TRUE)
p1 <- as.data.frame(t(p1))
p1 <- as.matrix(p1)
m1 <- metab %>% select(all_of(names(d4)))
m1 <- as.data.frame(t(m1))
m1 <- scale(m1, scale=TRUE, center=TRUE)
m1 <- as.data.frame(t(m1))
m1 <- as.matrix(m1)
aab1 <- aab %>% select(all_of(names(d4)))
aab1 <- as.data.frame(t(aab1))
aab1 <- scale(aab1, scale=TRUE, center=TRUE)
aab1 <- as.data.frame(t(aab1))
aab1 <- as.matrix(aab1)
c1 <- ciber %>% select(all_of(names(d4)))
c1 <- as.data.frame(t(c1))
c1 <- scale(c1, scale=TRUE, center=TRUE)
c1 <- as.data.frame(t(c1))
c1 <- as.matrix(c1)


#this has recuded the number to 468 for each
#make it into a list
corrected_omics_list2 <- list(t1, p1, m1, aab1, c1)
##################################################################################################
#do the MOFA
#get the code working for the initial MOFA
MOFAobject <- create_mofa(corrected_omics_list2)

#gives warning of duplicates features (probably proteomics and RNA-seq)
#should be okay as MOFA renames these
#NOTE: consider log-transforming/normalising the CIBER (albeit MOFA will scale too)
#CIBER was log-transformed for correction for residuals

#visualise MOFA
plot_data_overview(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)
head(model_opts)
#change below for different number of factors (normally 15)
#model_opts$num_factors <- 10
train_opts <- get_default_training_options(MOFAobject)
head(train_opts)
#default seed is 42 (no need to change for now)

data_opts <- get_default_data_options(MOFAobject)
head(data_opts)
#I do want to scale the views as everything is on a differing scale!
data_opts$scale_views <- TRUE

#now train the MOFA model!!
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

outfile = file.path(getwd(),"model.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)

#Warning message:
#  In .quality_control(object, verbose = verbose) :
#  Factor(s) 1, 2 are strongly correlated with the total number of expressed features for at least one of your omics. Such factors appear when there are differences in the total 'levels' between your samples, *sometimes* because of poor normalisation in the preprocessing steps.
#this error has occured before

#now process the MOFA
#get overview of data
plot_data_overview(MOFAobject.trained)

#add metadata
MOFAobject.trained@samples_metadata
sum(MOFAobject.trained@dimensions$N)

#also get the factor scores
factor_scores <- get_factors(MOFAobject.trained)
factor_scores_validation <- as.data.frame(factor_scores)
factor_scores_validation$samples <- rownames(factor_scores_validation)

views_names(MOFAobject.trained) <- c("Transcriptomics", "Proteomics", "Metabolomics", "Autoantibodies", "Leukocytes")
plot_variance_explained(MOFAobject.trained, x="view", y="factor")



varp <- plot_variance_explained(MOFAobject.trained, x="view", y="factor") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Generate the plot
varp <- plot_variance_explained(MOFAobject.trained, x = "view", y = "factor")
print(varp)

# Extract actual factor levels from the plot
factors_used <- ggplot_build(varp)$layout$panel_params[[1]]$y$get_labels()

# Reverse their order
varp <- varp +
  scale_y_discrete(limits = rev(factors_used)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



# Display the plot
varp


varp <- varp +
  scale_fill_gradientn(
    colours = c(
      "white",    # 0%
      "#FFF3B0",  # 2%
      "#FDB863",  # 5%
      "#E66101",  # 10%
      "#B2182B",  # 15%
      "#67001F"   # 20%
    ),
    values = scales::rescale(c(0, 2, 5, 10, 15, 20)),
    limits = c(0, 20),
    name = "Var. (%)"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 1,
      vjust = 0.5
    )
  )

print(varp)

#also get weights
weights <- get_weights(MOFAobject.trained, views = "all", factors = "all")
rna_weight <- as.data.frame(weights$Transcriptomics)
proteome_weight <- as.data.frame(weights$Proteomics)
metab_weight <- as.data.frame(weights$Metabolomics)

lapply(proteome_weight, hist)
lapply(metab_weight, hist)

#now investigate factors
plot_top_weights(MOFAobject.trained,
                 view = "Transcriptomics",
                 factor = 1,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Metabolomics",
                 factor = 1,
                 nfeatures = 25
)


plot_top_weights(MOFAobject.trained,
                 view = "Transcriptomics",
                 factor = 2,
                 nfeatures = 25
)


plot_top_weights(MOFAobject.trained,
                 view = "Proteomics",
                 factor = 3,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Metabolomics",
                 factor = 3,
                 nfeatures = 25
)


plot_top_weights(MOFAobject.trained,
                 view = "Proteomics",
                 factor = 4,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Metabolomics",
                 factor = 4,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Transcriptomics",
                 factor = 5,
                 nfeatures = 25
)


plot_top_weights(MOFAobject.trained,
                 view = "Transcriptomics",
                 factor = 6,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Metabolomics",
                 factor = 7,
                 nfeatures = 25)

plot_top_weights(MOFAobject.trained,
                 view = "Proteomics",
                 factor = 7,
                 nfeatures = 25)


plot_top_weights(MOFAobject.trained,
                 view = "Transcriptomics",
                 factor = 8,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Metabolomics",
                 factor = 8,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Proteomics",
                 factor = 8,
                 nfeatures = 25
)


plot_top_weights(MOFAobject.trained,
                 view = "Transcriptomics",
                 factor = 9,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Transcriptomics",
                 factor = 10,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Transcriptomics",
                 factor = 11,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Proteomics",
                 factor = 11,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Transcriptomics",
                 factor = 12,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Proteomics",
                 factor = 12,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Metabolomics",
                 factor = 12,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Transcriptomics",
                 factor = 13,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Proteomics",
                 factor = 13,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Transcriptomics",
                 factor = 14,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Proteomics",
                 factor = 14,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Transcriptomics",
                 factor = 15,
                 nfeatures = 25
)

plot_top_weights(MOFAobject.trained,
                 view = "Proteomics",
                 factor = 15,
                 nfeatures = 25
)

# RJ - why on earth did you call this fact9??

fact9 <- rna_weight %>% filter(abs(Factor10) >0.5)
fact9$gene <- rownames(fact9)

ggplot(fact9, aes(x = reorder(gene, -Factor10), y = Factor10)) +
  geom_bar(stat = 'identity', fill = 'purple') +
  xlab('Gene') +
  ylab('Weighting of gene') +
  ggtitle('Gene weight for each transcript’s contribution to factor 10 - validation analysis') +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


######################################################################################
#Similar style plot as in main for supplemental validation 

make_weight_plots <- function(df, factors, modality, n_top = 25) {
  
  df$X <- sub("_view_[0-9]+$", "", df$X)
  n_use <- min(n_top, nrow(df))
  
  plots <- lapply(factors, function(f) {
    
    top_idx <- order(abs(df[[f]]), decreasing = TRUE)[seq_len(n_use)]
    
    plot_df <- df[top_idx, ] %>%
      transmute(Feature = X, Weight = .data[[f]]) %>%
      mutate(
        Feature = reorder(Feature, Weight),
        Direction = ifelse(Weight > 0, "Positive", "Negative"),
        label_x = ifelse(Weight > 0, -0.02, 0.02),
        label_h = ifelse(Weight > 0, 1, 0)
      )
    
    ggplot(plot_df, aes(x = Weight, y = Feature, colour = Direction)) +
      geom_vline(xintercept = 0, linetype = "dotted", colour = "grey30") +
      geom_segment(aes(x = 0, xend = Weight, yend = Feature), linewidth = 0.6) +
      geom_point(size = 2.5) +
      geom_text(aes(x = label_x, label = Feature, hjust = label_h),
                colour = "black", size = 3) +
      scale_colour_manual(values = c(Positive = "#B2182B", Negative = "#2166AC"),
                          guide = "none") +
      scale_x_continuous(breaks = seq(-1, 1, 0.5)) +
      coord_cartesian(xlim = c(-1, 1)) +
      labs(title = paste0(sub("Factor", "Factor ", f), " - top ", n_use, " ", modality, " weights"),
           x = "Weight", y = NULL) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major.y = element_blank(),
        axis.text.y        = element_blank(),
        axis.line.x        = element_line(colour = "black", linewidth = 0.5),
        legend.position    = "none"
      )
  })
  names(plots) <- factors
  plots
}

# ---- 2. Prepare the validation RNA weights (genes are in rownames of rna_weight) ----
rna_val <- rna_weight %>% rownames_to_column("X")

# ---- 3. Factor 10 plot ----
# Option A: top 25 genes by absolute weight (matches the other panels)
p10_rna <- make_weight_plots(rna_val, "Factor10",
                             "transcriptomic (validation)")[["Factor10"]]

print(p10_rna)

# check that no weights fall outside the plotted x-axis range (-1, 1)
max(abs(rna_val$Factor10))

# ---- 4. Save ----
# height is fixed at 7; increase it if you show many genes (Option B)
ggsave("weights_Factor10_transcriptomic_validation.png", p10_rna,
       width = 6, height = 7, dpi = 300, bg = "white")

######################################################################################
#quickly check survival!
#get survival differences
#merge with MOFA datamort_df
#merge with MOFA data
mdf <- merge(mort_df, factor_scores_validation, by.x='id', by.y='samples')


sobj <- Surv(mdf$surv_time, mdf$event, type='right')
cox <- coxph(sobj ~ group1.Factor1 + group1.Factor2 + group1.Factor3 + group1.Factor4 + group1.Factor5 + group1.Factor6 + group1.Factor7 + group1.Factor8 + group1.Factor9 + group1.Factor10 + group1.Factor11 + group1.Factor12 + group1.Factor13+ group1.Factor14 + group1.Factor15, data=mdf)

ggforest(cox, data=mdf)

#now select those of interest
cox <- coxph(sobj ~ group1.Factor1 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + cbt_inflammation_crp_mgpl+ bs_bmi, data=mdf)
ggforest(cox, data=mdf)


cox <- coxph(sobj ~ group1.Factor2 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + cbt_inflammation_crp_mgpl+ bs_bmi, data=mdf)
ggforest(cox, data=mdf)

cox <- coxph(sobj ~  group1.Factor3 + age_diagnosis + sex+ hb_pvr_calc + hb_rap_m + cbt_inflammation_crp_mgpl+ bs_bmi, data=mdf)
ggforest(cox, data=mdf)

cox <- coxph(sobj ~ group1.Factor4 + age_diagnosis + sex+ hb_pvr_calc + hb_rap_m + cbt_inflammation_crp_mgpl+ bs_bmi, data=mdf)
ggforest(cox, data=mdf)

cox <- coxph(sobj ~ group1.Factor5 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + cbt_inflammation_crp_mgpl+ bs_bmi, data=mdf)
ggforest(cox, data=mdf)

cox <- coxph(sobj ~ group1.Factor6 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + cbt_inflammation_crp_mgpl+ bs_bmi, data=mdf)
ggforest(cox, data=mdf)

cox <- coxph(sobj ~ group1.Factor7 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + cbt_inflammation_crp_mgpl+ bs_bmi, data=mdf)
ggforest(cox, data=mdf)

cox <- coxph(sobj ~ group1.Factor8 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + cbt_inflammation_crp_mgpl+ bs_bmi, data=mdf)
ggforest(cox, data=mdf)

cox <- coxph(sobj ~ group1.Factor9 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + cbt_inflammation_crp_mgpl + bs_bmi, data=mdf)
ggforest(cox, data=mdf)

cox <- coxph(sobj ~ group1.Factor10 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + cbt_inflammation_crp_mgpl + bs_bmi + time_to_sampling, data=mdf)
summary(cox)
ggforest(cox, data=mdf)

cox <- coxph(sobj ~ group1.Factor10 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + cbt_inflammation_crp_mgpl + bs_bmi, data=mdf)
summary(cox)
ggforest(cox, data=mdf)

###############################################
#visually nice CoxPH

# Create a copy of the data with nicer column names for plotting
mdf_labels <- mdf
colnames(mdf_labels)[colnames(mdf_labels) == "group1.Factor10"] <- "Factor 10"
colnames(mdf_labels)[colnames(mdf_labels) == "age_diagnosis"] <- "Age at diagnosis"
colnames(mdf_labels)[colnames(mdf_labels) == "sex"] <- "Sex"
colnames(mdf_labels)[colnames(mdf_labels) == "hb_pvr_calc"] <- "PVR"
colnames(mdf_labels)[colnames(mdf_labels) == "hb_rap_m"] <- "RAP"
colnames(mdf_labels)[colnames(mdf_labels) == "cbt_inflammation_crp_mgpl"] <- "CRP"
colnames(mdf_labels)[colnames(mdf_labels) == "bs_bmi"] <- "BMI"

# Refit the same model on the renamed data (just for plotting)
cox_labels <- coxph(sobj ~ `Factor 10` + `Age at diagnosis` + Sex + PVR + RAP + CRP + BMI, data = mdf_labels)

# Plot with nice labels
ggforest(cox_labels, data = mdf_labels)

#######################################################################################


cox <- coxph(sobj ~group1.Factor11 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + cbt_inflammation_crp_mgpl + bs_bmi, data=mdf)
ggforest(cox, data=mdf)

cox <- coxph(sobj ~ group1.Factor12 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + cbt_inflammation_crp_mgpl + bs_bmi, data=mdf)
ggforest(cox, data=mdf)

cox <- coxph(sobj ~ group1.Factor13 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + cbt_inflammation_crp_mgpl + bs_bmi, data=mdf)
ggforest(cox, data=mdf)

cox <- coxph(sobj ~ group1.Factor14 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + cbt_inflammation_crp_mgpl + bs_bmi, data=mdf)
ggforest(cox, data=mdf)

cox <- coxph(sobj ~ group1.Factor15 + age_diagnosis + sex + hb_pvr_calc + hb_rap_m + cbt_inflammation_crp_mgpl + bs_bmi, data=mdf)
ggforest(cox, data=mdf)


#get correlation plot

#correlate with original factors
#definitive_mofa_model_factor_scores_adjused_for_filtering <- readRDS("/mnt/session_data/definitive_mofa_model_factor_scores_adjused_for_filtering.rds")
definitive_mofa_model_May25_PC_WBC_resids <- readRDS("/mnt/session_data/definitive_mofa_model_May25_PC_WBC_resids.rds")
definitive_mofa_model_factor_scores_adjused_for_filtering <- definitive_mofa_model_May25_PC_WBC_resids
#get factor correlations
#now get correlations between factors
colnames(definitive_mofa_model_factor_scores_adjused_for_filtering) <- ifelse(colnames(definitive_mofa_model_factor_scores_adjused_for_filtering) == "samples", "samples", paste0("definitive_", colnames(definitive_mofa_model_factor_scores_adjused_for_filtering)))
colnames(factor_scores_validation) <- ifelse(colnames(factor_scores_validation) == "samples", "samples", paste0("compl_w_miss__", colnames(factor_scores_validation)))
#now names are changes, plot in a df for correlations
df1 <- merge(definitive_mofa_model_factor_scores_adjused_for_filtering, factor_scores_validation, by='samples')
df1 <- df1[,-1]

#now plot!
row_vars <- grep("^definitive_", colnames(df1), value = TRUE)
col_vars <- setdiff(colnames(df1), row_vars)

# Subset the dataframe
df_rows <- df1[, row_vars]
df_cols <- df1[, col_vars]

library(seriation)

#Compute correlation matrix
cor_matrix <- cor(df1[, row_vars], df1[, col_vars], method = "pearson")
abs_cor <- abs(cor_matrix)

#For each row, find the best matching column (max correlation)
best_matches <- apply(abs_cor, 1, which.max)

# Sort rows by the index of their best match
row_order <- order(best_matches)
col_order <- best_matches[row_order]  

# Reorder correlation matrix accordingly
cor_matrix_reordered <- cor_matrix[row_order, col_order]

# Plot
pheatmap(cor_matrix_reordered,cluster_rows = FALSE,cluster_cols = FALSE,display_numbers = TRUE)

##############################################################
#cleaner plot
# Clean up row names (definitive_group1.FactorX -> Factor X)
rownames(cor_matrix_reordered) <- gsub(".*Factor", "Factor ", rownames(cor_matrix_reordered))
# Clean up column names (compl_w_miss__group1.FactorX -> Validation Factor X)
colnames(cor_matrix_reordered) <- gsub(".*Factor", "Validation Factor ", colnames(cor_matrix_reordered))
# Plot
pheatmap(cor_matrix_reordered, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE)
#===============================================================================================================================
#associate factors scores of first itteration of MOFA with REVEAL score and with incident/prevalent status
incident_and_prevalent_cases <- readRDS("/mnt/mounted-data/incident_and_prevalent_cases.rds")
REVEAL_LITE_score <- readRDS("/mnt/mounted-data/REVEAL_LITE_score.rds")


reveal <- merge(REVEAL_LITE_score, definitive_mofa_model_factor_scores_adjused_for_filtering, by.x='id', by.y='samples')
rev_keep <- reveal[,1:2]
#157/179 have enough data for REVEAL-lite

#whilst mostly normal, not entirely --> do both spearman and pearson...
cor_results <- data.frame(variable = character(), correlation = numeric(), p_value = numeric())
reveal <- reveal[,-1]

for (var in names(reveal)[!names(reveal) %in% "reveal_lite"]) {
  test <- cor.test(reveal[[var]], reveal$reveal_lite, method = "spearman")
  cor_results <- rbind(cor_results, 
                       data.frame(variable = var, correlation = test$estimate, p_value = test$p.value))
}

# Prepare matrix for pheatmap
cor_matrix <- matrix(cor_results$correlation, nrow = 1, 
                     dimnames = list("reveal_lite", cor_results$variable))
pval_matrix <- matrix(formatC(cor_results$p_value, 4, digits = 2), nrow = 1, 
                      dimnames = list("reveal_lite", cor_results$variable))

# Customize the annotation of p-values in the cells
cell_text <- matrix(paste0("p=", pval_matrix), nrow = 1, 
                    dimnames = list("reveal_lite", cor_results$variable))

# Rename the column of cor_matrix to "REVEAL lite"
rownames(cor_matrix) <- "REVEAL lite"

colnames(cor_matrix) <- c("Factor 1",  "Factor 2",  "Factor 3",  "Factor 4",  "Factor 5", "Factor 6" , "Factor 7" , "Factor 8",  "Factor 9" , "Factor 10", "Factor 11", "Factor 12", "Factor 13", "Factor 14", "Factor 15")


# Plot the heatmap
pheatmap(cor_matrix, 
         display_numbers = cell_text,
         number_color = "black", # Color of p-values
         color = colorRampPalette(c("blue", "white", "red"))(50), # Color scale
         fontsize_number = 10, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         breaks = seq(-1, 1, length.out = 51), 
         main ='Spearman correlation of factor scores with REVEAL lite',
         legend_breaks = c(-1, -0.5, 0, 0.5, 1))

############################################################################################
#now test with incident/prevalent proportions
#Given roughly normal distribution, ANOVA should do it for now
#decide on this later though
ip <- merge(incident_and_prevalent_cases, definitive_mofa_model_factor_scores_adjused_for_filtering, by.x='id', by.y='samples')
library(Publish)
library(Hmisc)
#use 12 month for safey
u1 <- univariateTable(prevalent12months ~ definitive_group1.Factor1 + definitive_group1.Factor2 + definitive_group1.Factor3 + definitive_group1.Factor4 + definitive_group1.Factor5 + definitive_group1.Factor6 + definitive_group1.Factor7 + definitive_group1.Factor8 + definitive_group1.Factor9 + definitive_group1.Factor10 + definitive_group1.Factor11 + definitive_group1.Factor12 + definitive_group1.Factor13 + definitive_group1.Factor14 + definitive_group1.Factor15, data=ip, column.percent = TRUE, show.totals=TRUE)
u2 <- summary(u1)

#also check vasoresponders! 
vaso <- clin_dat_at_diagnosis_clean %>% select(id, hv_vasodilator_responder)
vaso$hv_vasodilator_responder <- as.character(vaso$hv_vasodilator_responder)
vaso$hv_vasodilator_responder <- ifelse(is.na(vaso$hv_vasodilator_responder), 'no', vaso$hv_vasodilator_responder)

names(vaso) <- c('id', 'vasoresponder')
vaso$vasoresponder <- as.factor(vaso$vasoresponder)


ip2 <- merge(vaso, definitive_mofa_model_factor_scores_adjused_for_filtering, by.x='id', by.y='samples')
u1 <- univariateTable(vasoresponder ~ definitive_group1.Factor1 + definitive_group1.Factor2 + definitive_group1.Factor3 + definitive_group1.Factor4 + definitive_group1.Factor5 + definitive_group1.Factor6 + definitive_group1.Factor7 + definitive_group1.Factor8 + definitive_group1.Factor9 + definitive_group1.Factor10 + definitive_group1.Factor11 + definitive_group1.Factor12 + definitive_group1.Factor13 + definitive_group1.Factor14 + definitive_group1.Factor15, data=ip2, column.percent = TRUE, show.totals=TRUE)
u2 <- summary(u1)

clin <- merge(clin, rev_keep, by='id', all=T)
clin <- merge(clin, time_df[,c(1,3)], by='id', all=T)
clin$reveal_lite <- clin$reveal_lite.x
clin <- clin %>% select(!any_of(c('reveal_lite.x', 'reveal_lite.y')))
clin$time_to_sampling <- clin$time_to_sampling.x
clin <- clin %>% select(!any_of(c('time_to_sampling.x', 'time_to_sampling.y')))
#===========================================
#first with small dataset and HC clustering
clin_small <- clin %>% select(id, hb_pvr_calc, hb_pap_m, hb_pawp_m, hb_rap_m, cbt_card_ntprobnp_ngpl, reveal_lite)

#merge with MOFA
merged_df <- merge(clin_small, definitive_mofa_model_factor_scores_adjused_for_filtering, by.x = 'id', by.y = 'samples')

num_cols_x <- names(clin_small)[sapply(clin_small, is.numeric)]
num_cols_y <- names(definitive_mofa_model_factor_scores_adjused_for_filtering)[sapply(definitive_mofa_model_factor_scores_adjused_for_filtering, is.numeric)]

#do spearman
cor_results <- rcorr(as.matrix(merged_df[, c(num_cols_x, num_cols_y)]), type = "spearman")

cor_matrix <- cor_results$r[num_cols_x, num_cols_y]
p_matrix <- cor_results$P[num_cols_x, num_cols_y]

#pval formatting
formatted_p_values <- ifelse(p_matrix < 0.001, "***",
                             ifelse(p_matrix < 0.01, "**",
                                    ifelse(p_matrix < 0.05, "*", "")))

cor_df <- as.data.frame(as.table(cor_matrix))
pval_df <- as.data.frame(as.table(p_matrix))
cor_df$pval <- pval_df$Freq
cor_df$label <- paste0(round(cor_df$Freq, 2), formatted_p_values)

#add HC (no dendrograms and same plot style for consistency)
row_order <- hclust(dist(cor_matrix))$order
col_order <- hclust(dist(t(cor_matrix)))$order

cor_df$Var1 <- factor(cor_df$Var1, levels = rownames(cor_matrix)[row_order])
cor_df$Var2 <- factor(cor_df$Var2, levels = colnames(cor_matrix)[col_order])

ggplot(cor_df, aes(Var2, Var1, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = label), color = "black", size = 3) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-1, 1)
  ) +
  theme_pubr() +
  labs(x = "MOFA Factors", y = "Clinical Variables", fill = "Spearman R") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


display_mat <- matrix(
  paste0(round(cor_matrix, 2), formatted_p_values), 
  nrow = nrow(cor_matrix)
)

# 2. Plot
pheatmap(cor_matrix, 
         display_numbers = display_mat, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         clustering_distance_rows = "euclidean", # matches your hclust(dist())
         clustering_distance_cols = "euclidean",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-1, 1, length.out = 101),
         main = "MOFA Factor correlations with clinical parameters",
         fontsize_number = 8)

######################################################################
#Make matrix visually clearer:
# --- Clean up axis labels ---

# Y-axis: map clinical variable codes to readable names
clin_label_map <- c(
  cbt_card_ntprobnp_ngpl = "NT-proBNP",
  reveal_lite             = "REVEAL Lite 2",
  hb_pvr_calc              = "PVR",
  hb_pap_m                 = "mPAP",
  hb_pawp_m                = "PAWP",
  hb_rap_m                 = "mRAP"
)

row_labels <- clin_label_map[rownames(cor_matrix)]

# X-axis: strip "definitive_group1." prefix and add space before number
col_labels <- gsub("^definitive_group1\\.", "", colnames(cor_matrix))
col_labels <- gsub("Factor", "Factor ", col_labels)

# 2. Plot
pheatmap(cor_matrix, 
         display_numbers = display_mat, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-1, 1, length.out = 101),
         main = "MOFA Factor correlations with clinical parameters",
         fontsize_number = 8,
         labels_row = row_labels,
         labels_col = col_labels)

#===================================================================
#now test associations with a variety of clinical points
merged_df <- merge(clin, definitive_mofa_model_factor_scores_adjused_for_filtering, 
                   by.x = 'id', by.y = 'samples')

num_cols_x <- names(clin)[sapply(clin, is.numeric)]
num_cols_y <- names(definitive_mofa_model_factor_scores_adjused_for_filtering)[sapply(definitive_mofa_model_factor_scores_adjused_for_filtering, is.numeric)]



cor_results <- rcorr(as.matrix(merged_df[, c(num_cols_x, num_cols_y)]), type = "spearman")

cor_matrix <- cor_results$r[num_cols_x, num_cols_y]


cor_matrix <- cor_results$r[num_cols_x, num_cols_y]
p_matrix   <- cor_results$P[num_cols_x, num_cols_y]

p_vec     <- as.vector(p_matrix)                      
fdr_vec   <- p.adjust(p_vec, method = "fdr")         
fdr_matrix <- matrix(fdr_vec,
                     nrow = nrow(p_matrix),
                     ncol = ncol(p_matrix),
                     dimnames = dimnames(p_matrix))    

cor_results$fdr <- fdr_matrix

formatted_p_values <- ifelse(fdr_matrix < 0.001, "***",
                             ifelse(fdr_matrix < 0.01, "**",
                                    ifelse(fdr_matrix < 0.05, "*", "")))

cor_df <- as.data.frame(as.table(cor_matrix))
fdr_df <- as.data.frame(as.table(fdr_matrix))
cor_df$fdr <- fdr_df$Freq
cor_df$label <- paste0(round(cor_df$Freq, 2), formatted_p_values)

#now plot
ggplot(cor_df, aes(Var2, Var1, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = label), color = "black") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-1, 1)
  ) +
  theme_pubr() +
  labs(x = "MOFA Factors", y = "Clinical Variables", fill = "Spearman R") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#####################################################################
#clean plot:

#Y-axis: map clinical variable codes to readable names
#NOTE: double-check these names against rownames(cor_matrix) / num_cols_x
clin_label_map <- c(
  time_to_sampling              = "Time to Sampling",
  reveal_lite                   = "REVEAL Lite 2",
  cbt_thyr_tsh_mupl              = "TSH",
  lf_kco_pc                     = "KCO",
  lf_fvc_pc                     = "FVC",
  lf_fev1_pc                    = "FEV1",
  ep_1_distance_meters          = "6MWD",
  hb_cardiac_index_value_1      = "Cardiac Index",
  hb_cardiac_output_value_1     = "Cardiac Output",
  hb_pawp_m                     = "PAWP",
  hb_pap_m                      = "mPAP",
  hb_rap_m                      = "mRAP",
  hb_pvr_calc                    = "PVR",
  cbt_haem_hb_gpl                = "Haemoglobin",
  cbt_haem_wbc_x10e9pl            = "WBC",
  cbt_haem_redcell_dist_pct       = "RDW",
  cbt_haem_platelets_x10e9pl      = "Platelets",
  cbt_inflammation_crp_mgpl       = "CRP",
  cbt_card_ntprobnp_ngpl          = "NT-proBNP",
  bs_bmi                        = "BMI",
  age_diagnosis                 = "Age at Diagnosis"
)

# X-axis: strip "group1." prefix and add space before number
col_label_map <- setNames(
  gsub("Factor", "Factor ", gsub("^group1\\.", "", colnames(cor_matrix))),
  colnames(cor_matrix)
)

# --- Plot with cleaned labels ---
ggplot(cor_df, aes(Var2, Var1, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = label), color = "black") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-1, 1)
  ) +
  scale_x_discrete(labels = col_label_map) +
  scale_y_discrete(labels = clin_label_map) +
  theme_pubr() +
  labs(x = "MOFA Factors", y = "Clinical Variables", fill = "Spearman R") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))



# --- Add (n=) of non-missing values to y-axis labels ---

# Non-missing count per clinical variable (in the merged dataset)
n_vec <- colSums(!is.na(merged_df[, num_cols_x]))

# Readable name if mapped, otherwise fall back to the raw variable name
base_names <- ifelse(is.na(clin_label_map[names(n_vec)]),
                     names(n_vec),
                     clin_label_map[names(n_vec)])

clin_label_map_n <- setNames(paste0(base_names, " (n=", n_vec, ")"),
                             names(n_vec))

# --- Plot ---
ggplot(cor_df, aes(Var2, Var1, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = label), color = "black") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-1, 1)
  ) +
  scale_x_discrete(labels = col_label_map) +
  scale_y_discrete(labels = clin_label_map_n) +
  theme_pubr() +
  labs(x = "MOFA Factors", y = "Clinical Variables", fill = "Spearman R") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


# X-axis: strip "definitive_group1." prefix and add space before number
col_label_map <- setNames(
  gsub("Factor", "Factor ", gsub("^definitive_group1\\.", "", colnames(cor_matrix))),
  colnames(cor_matrix)
)

# --- Plot ---
ggplot(cor_df, aes(Var2, Var1, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = label), color = "black") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-1, 1)
  ) +
  scale_x_discrete(labels = col_label_map) +
  scale_y_discrete(labels = clin_label_map_n) +
  theme_pubr() +
  labs(x = "MOFA Factors", y = "Clinical Variables", fill = "Spearman R") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

##########################################################################
#Make table
cor_table <- data.frame(Variable1 = rep(rownames(cor_matrix), ncol(cor_matrix)),
                        Variable2 = rep(colnames(cor_matrix), each = nrow(cor_matrix)),
                        Spearman_R = as.vector(cor_matrix),
                        P_value = as.vector(p_matrix))

cor2 <- cor_table
cor2$padj <- p.adjust(cor2$P_value, method='fdr')
cor2 <- cor2 %>% filter(padj <0.05)



#repeat difference for sex and disease type
ip3 <- merge(clin_dat_at_diagnosis_clean, definitive_mofa_model_May25_PC_WBC_resids, by.x='id', by.y='samples')
u1 <- univariateTable(sex ~ group1.Factor1 + group1.Factor2 + group1.Factor3 + group1.Factor4 + group1.Factor5 + group1.Factor6 + group1.Factor7 + group1.Factor8 + group1.Factor9 + group1.Factor10 + group1.Factor11 + group1.Factor12 + group1.Factor13 + group1.Factor14 + group1.Factor15, data=ip3, column.percent = TRUE, show.totals=TRUE)
u2 <- summary(u1)

ip3$diagnosis_verified <- droplevels(ip3$diagnosis_verified)
u1 <- univariateTable(diagnosis_verified ~ group1.Factor1 + group1.Factor2 + group1.Factor3 + group1.Factor4 + group1.Factor5 + group1.Factor6 + group1.Factor7 + group1.Factor8 + group1.Factor9 + group1.Factor10 + group1.Factor11 + group1.Factor12 + group1.Factor13 + group1.Factor14 + group1.Factor15, data=ip3, column.percent = TRUE, show.totals=TRUE)
u2 <- summary(u1)



#now visualise
ip4 <- ip3 %>% pivot_longer(304:318, names_to = 'Factor', values_to = 'factor_value')

ggplot(ip4, aes(x=factor_value, y=diagnosis_verified, fill = diagnosis_verified)) +
  geom_boxplot() +
  facet_wrap(~Factor) +
  ggtitle('Factor scores per factor for each diagnosis') +
  xlab('Factor value') +
  ylab('Diagnosis') +
  theme_pubr()


###########################################################################################3
#now correlate MOFA factors with autoabs
#definitive_mofa_model_factor_scores_adjused_for_filtering <- readRDS("/mnt/session_data/definitive_mofa_model_factor_scores_adjused_for_filtering.rds")
definitive_mofa_model_May25_PC_WBC_resids <- readRDS("/mnt/session_data/definitive_mofa_model_May25_PC_WBC_resids.rds")
definitive_mofa_model_factor_scores_adjused_for_filtering <- definitive_mofa_model_May25_PC_WBC_resids

Autoantibody_levels_patients_log_transformed <- readRDS("/mnt/session_data/mounted-data-readonly/PAH-omics/omics/Autoantibodies/Autoantibody_levels_patients_log_transformed.rds")
follow_up_data <- readRDS("/mnt/session_data/mounted-data-readonly/PAH-OpenClinica/20240726-freeze/follow_up_data.rds")

#get labels right
fu <- follow_up_data %>% select(id, id_cohort) %>% unique()

Autoantibody_levels_patients_log_transformed <- merge(Autoantibody_levels_patients_log_transformed, fu, by.x='sample', by.y='id_cohort')
aabs <- merge(Autoantibody_levels_patients_log_transformed, definitive_mofa_model_factor_scores_adjused_for_filtering, by.x='id', by.y='samples')

#now get Spearman cor 
num_cols_x <- names(Autoantibody_levels_patients_log_transformed)[sapply(Autoantibody_levels_patients_log_transformed, is.numeric)]
num_cols_y <- names(definitive_mofa_model_factor_scores_adjused_for_filtering)[sapply(definitive_mofa_model_factor_scores_adjused_for_filtering, is.numeric)]

library(Hmisc)

cor_results <- rcorr(as.matrix(aabs[, c(num_cols_x, num_cols_y)]), type = "spearman")


cor_matrix <- cor_results$r[num_cols_x, num_cols_y]
p_matrix <- cor_results$P[num_cols_x, num_cols_y]


formatted_p_values <- ifelse(p_matrix < 0.001, "***", 
                             ifelse(p_matrix < 0.01, "**", 
                                    ifelse(p_matrix < 0.05, "*", "")))

cor_df <- as.data.frame(as.table(cor_matrix))
pval_df <- as.data.frame(as.table(p_matrix))


cor_df$pval <- pval_df$Freq
cor_df$label <- paste0(round(cor_df$Freq, 2), formatted_p_values)

#now plot
ggplot(cor_df, aes(Var2, Var1, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = label), color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_minimal() +
  labs(x = "MOFA Factors", y = "Autoantibodies", fill = "Spearman R")

#Make table
cor_table <- data.frame(Variable1 = rep(rownames(cor_matrix), ncol(cor_matrix)),
                        Variable2 = rep(colnames(cor_matrix), each = nrow(cor_matrix)),
                        Spearman_R = as.vector(cor_matrix),
                        P_value = as.vector(p_matrix))

cor2 <- cor_table
cor2$padj <- p.adjust(cor2$P_value, method='fdr')
cor2 <- cor2 %>% filter(padj <0.05)

#=========================================================================================================
#EdB continued here on 12/02/25
#now correlate all factor scores in the new df
facts <- definitive_mofa_model_factor_scores_adjused_for_filtering
facts <- facts %>% select(-samples)

MOFA_trascriptomic_weights_6Nov_25 <- read.csv2("/mnt/session_data/MOFA_trascriptomic_weights_6Nov_25.csv")

#get Spearman
cor_matrix <- cor(facts, method = "pearson", use = "pairwise.complete.obs")

color_palette <- colorRampPalette(c("blue", "white", "red"))(99)
breaks_list <- seq(-1, 1, length.out = 100)  

pheatmap(cor_matrix, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, breaks=breaks_list, color=color_palette)  

plot(definitive_mofa_model_factor_scores_adjused_for_filtering$group1.Factor10, definitive_mofa_model_factor_scores_adjused_for_filtering$group1.Factor12)

#also plot the weights against each other
plot(rna_weight$Factor10, rna_weight$Factor12)


rna3 <- MOFA_trascriptomic_weights_6Nov_25 %>% filter(abs(Factor11) > 0.25 | abs(Factor13) > 0.25)
rownames(rna3) <- rna3$X

plot(rna3$Factor11, rna3$Factor13)

#now annotate
cor.test(MOFA_trascriptomic_weights_6Nov_25$Factor11, MOFA_trascriptomic_weights_6Nov_25$Factor13, method='spearman')

plot(rna3$Factor11, rna3$Factor13, 
     xlab = "Factor 11", 
     ylab = "Factor 13", 
     main = "Scatter Plot of Factor 11 vs Factor 13", 
     pch = 16, col = "green")

text(rna3$Factor11, rna3$Factor13, labels = rownames(rna3), pos = 4, cex = 0.7, col = "blue")


#EdB continued on 13/02/25
#check green genes from Chris

rna4 <- rna_weight %>% filter(Factor13 < 0 & Factor12 < 0.2)


plot(rna4$Factor11, rna4$Factor13, 
     xlab = "Factor 11", 
     ylab = "Factor 13", 
     main = "Scatter Plot of Factor 11 vs Factor 13", 
     pch = 16, col = "black")

text(rna4$Factor10, rna4$Factor12, labels = rownames(rna4), pos = 4, cex = 0.7, col = "red")

#now do genetic correlations
rna4$gene <- rownames(rna4)
#drop view 1
rna4[] <- lapply(rna4, function(x) if(is.character(x)) sub("_view_1", "", x) else x)

#extract this from the RNA data
rna5 <- transcript %>% filter(rownames(transcript) %in% rna4$gene)

#remove cols with NA
rna5 <- rna5 %>% select(where(~ sum(is.na(.)) < nrow(rna5)))

#get correlations
cor_matrix <- cor(as.data.frame(t(rna5)), method = "pearson", use = "pairwise.complete.obs")

#Pheatmap
pheatmap(cor_matrix, 
         main = "Pearson Correlation Heatmap", 
         cluster_rows = TRUE,  
         cluster_cols = TRUE,  
         clustering_method = "complete",  
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         border_color = NA,
         breaks = seq(-1, 1, length.out = 101))

#now also get the others
rna6 <- rna_weight %>% filter(Factor10 > 0.3 & Factor12 < -0.3)

rna6$gene <- rownames(rna6)
#drop view 1
rna6[] <- lapply(rna6, function(x) if(is.character(x)) sub("_view_1", "", x) else x)

#extract this from the RNA data
rna6 <- transcript %>% filter(rownames(transcript) %in% rna6$gene)

rna6 <- rna6 %>% select(where(~ sum(is.na(.)) < nrow(rna6)))


#get correlations
# Assign row names
rownames(rna5) <- paste0("green_genes_", 1:nrow(rna5))
rownames(rna6) <- paste0("blue_genes_", 1:nrow(rna6))

# Combine both datasets
combined <- rbind(rna5, rna6)

# Compute correlation matrix
cor_matrix <- cor(t(combined))

# Fix annotation row creation
annotation_row <- data.frame(Source = c(rep("green", nrow(rna5)), rep("blue", nrow(rna6))))
rownames(annotation_row) <- rownames(combined)  # Now it matches the combined rownames

ann_colors <- list(Source = c(green = "green", blue = "blue"))

# Generate heatmap
pheatmap(cor_matrix, 
         annotation_row = annotation_row, 
         annotation_colors = ann_colors, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         show_rownames = FALSE, 
         show_colnames = FALSE,
         breaks = seq(-1, 1, length.out = 100),
         main = "Correlation Heatmap")


#####################################################################################################################
#full DGE for IFN genes
#set wd
setwd('Myfilepath/Interferone response')

#load the packages
library(tidyverse)
library(edgeR)
library(broom)
library("EnhancedVolcano")
library(readxl)
library(ggpubr)
library(pheatmap)
library(factoextra)
library(FactoMineR)

#load in the data
interfereone_repsonse_genes_rna_seq <- readRDS("Myfilepath/Interferone response/interfereone_repsonse_genes_rna_seq.rds")
genes_keep <- colnames(interfereone_repsonse_genes_rna_seq[,5:225])
write.csv2(genes_keep, file='list_of_IFN_response.csv')
PCs_and_WBCs_RNA_seq <- readRDS("Myfilepath/Omics/PCs_and_WBCs_RNA_seq.rds")
Complete_RNA_seq_cohort_pheno <- readRDS("Myfilepath/Omics/Complete_RNA_seq_cohort_pheno.rds")
Complete_RNA_seq_cohort_counts <- readRDS("Myfilepath/Omics/Complete_RNA_seq_cohort_counts.rds")


#prepare the data (incl.removal of duplicates)
Complete_RNA_seq_cohort_pheno$SampleID <- as.character(Complete_RNA_seq_cohort_pheno$SampleID)
#Complete_RNA_seq_cohort_pheno$SampleID <- ifelse(Complete_RNA_seq_cohort_pheno$SampleID == OMITTED HERE FOR PRIVACY', Complete_RNA_seq_cohort_pheno$SampleID)

doubles <- Complete_RNA_seq_cohort_pheno %>% filter(Repeats == 2) %>% group_by(External.ID) %>% filter(n() > 1) %>% slice(2) %>%  ungroup()
counts <- Complete_RNA_seq_cohort_counts %>% select(!all_of(doubles$SampleID))
counts_removed <- Complete_RNA_seq_cohort_counts %>% select(all_of(doubles$SampleID))

#this makes sense, duplicates removed and number of samples also correct. 
sampleinfo <- Complete_RNA_seq_cohort_pheno %>% filter(SampleID %in% colnames(counts))
sampleinfo <- sampleinfo %>% select(SampleID, diagnosis_verified, PAH, Description.1, External.ID, CohortID, OpenClinicaID) %>% unique()
#now missing one, which has the added _
#sampleinfo[504,1] <- OMITTED FOR PRIVACXY
sampleinfo[504,2] <- 'PAH'

#remove all the relatives
sampleinfo <- sampleinfo %>% filter(Description.1 !='RELATIVE_1')
summary(as.factor(sampleinfo$diagnosis_verified))

sampleinfo$PAH <- ifelse(is.na(sampleinfo$PAH) & sampleinfo$diagnosis_verified == '0', 0, sampleinfo$PAH)
sampleinfo$PAH <- ifelse(is.na(sampleinfo$PAH), 1, sampleinfo$PAH)
sampleinfo <- unique(sampleinfo)

sampleinfo <- sampleinfo %>% filter(diagnosis_verified != 'PVOD/PCH')

rownames(sampleinfo) <- sampleinfo$SampleID
summary(as.factor(sampleinfo$PAH))
# 0   1 
#114 376 

#bin CTEPH patient here
#sampleinfo <- sampleinfo %>% filter(SampleID != 'OMITTED HERE FOR PRIVACY')
#remove Sheffield PAH here as these are PAH, despite being recorded as controls
sampleinfo <- sampleinfo %>% filter(Description.1 != 'Shef PAH')
summary(as.factor(sampleinfo$PAH))

rna_req_to_keep <- sampleinfo[,c(1,5,6,7,3)]
rna_seq_to_keep <- rna_req_to_keep
rna_seq_to_keep$OpenClinicaID <- as.character(rna_seq_to_keep$OpenClinicaID)
rna_seq_to_keep$OpenClinicaID <- ifelse(rna_seq_to_keep$CohortID == 'OMITTED FOR PRIVACY', 'OMITTED FOR PRIVACY', rna_seq_to_keep$OpenClinicaID)
rna_seq_to_keep$OpenClinicaID <- ifelse(rna_seq_to_keep$CohortID == 'OMITTED FOR PRIVACY', 'OMITTED FOR PRIVACY', rna_seq_to_keep$OpenClinicaID)


control_basics <- read_excel("Myfilepath/Omics/control basics.xlsx")
pheno_data_Allan_Mark_longitudinal <- read_excel("Myfilepath/Omics/pheno_data_Allan_Mark_longitudinal.xls")

phen1 <- pheno_data_Allan_Mark_longitudinal %>% select(SampleID, ethnicity, Age.atSample, Age_controls_sample_patients_diagnosis, sex)
phen1$Age_controls_sample_patients_diagnosis <- ifelse(is.na(phen1$Age_controls_sample_patients_diagnosis), phen1$Age.atSample, phen1$Age_controls_sample_patients_diagnosis)

phen1 <- phen1 %>% select(SampleID, Age_controls_sample_patients_diagnosis, sex, ethnicity)
names(phen1) <- c('SampleID', 'Age_at_sample', 'sex', 'ethnicity')
phen2 <- control_basics[,c(1,3,4,5)]
names(phen2) <- c('SampleID', 'Age_at_sample', 'sex', 'ethnicity')

phenodat <- unique(rbind(phen1, phen2))
phenodat[phenodat == "NA"] <- NA

#now ensure it only has those within the sample data
phenodat <- phenodat %>% filter(SampleID %in% rna_seq_to_keep$SampleID)
phenodat <- merge(rna_seq_to_keep, phenodat, by='SampleID', all=T)

#get control and patient data
clin_dat_at_diagnosis_clean <- readRDS("Myfilepath/Dataset for CYNAPSE/clin_dat_at_diagnosis_clean.rds")
follow_up_data <- readRDS("Myfilepath/Dataset for CYNAPSE/follow_up_data.rds")

clin_fu <- clin_dat_at_diagnosis_clean %>% select(id, sex, diagnosis_date, diagnosis_verified, age_diagnosis)
fu <- follow_up_data %>% select(id, visit, visit_date) %>% filter(visit ==1) %>% unique()

clin_fu <- merge(clin_fu, fu, by='id', all=T)

load("Myfilepath/2023-12-02/cleaned/data_clean.RData")

crtl_in <- data.clean.unrelctrls %>% filter(id %in% phenodat$OpenClinicaID)
relat_in <- data.clean.relative %>% filter(id %in% phenodat$OpenClinicaID)
relat_in <- relat_in %>% filter(visit == 1)
#relative is actually patient and phenoconverter, so only control FU and clin_fu needed

phenodat_miss <- phenodat[!complete.cases(phenodat), ]
crtl_in <- crtl_in %>% filter(id %in% phenodat_miss$OpenClinicaID)
clin_fu <- clin_fu %>% filter(id %in% phenodat_miss$OpenClinicaID)
clin_fu$timediff <- as.numeric(clin_fu$visit_date - clin_fu$diagnosis_date)/365
clin_fu$sampledate2 <- clin_fu$age_diagnosis + clin_fu$timediff
clin_fu <- clin_fu %>% select(id, sex, sampledate2)
names(clin_fu) <- c('id', 'sex2', 'sampledate2')

#repeat prep for controls
crtl_in$sampledate2 <- as.numeric(crtl_in$pah_control_visit_year) - as.numeric(crtl_in$pah_control_birth_year)
crtl_in <- crtl_in %>% select(id, pah_control_sex, sampledate2,pah_control_ethnic_category)

clin_fu$ethnicity2 <- NA

names(crtl_in) <- c('id', 'sex2', 'sampledate2', 'ethnicity2')
data_compl <- rbind(clin_fu, crtl_in)

phenodat <- merge(phenodat, data_compl, by.x='OpenClinicaID', by.y='id', all=T)
phenodat$Age_at_sample <- ifelse(is.na(phenodat$Age_at_sample), phenodat$sampledate2, phenodat$Age_at_sample)
phenodat$ethnicity <- ifelse(is.na(phenodat$ethnicity), phenodat$ethnicity2, phenodat$ethnicity)
phenodat$sex <- ifelse(is.na(phenodat$sex), phenodat$sex2, phenodat$sex)

phenodat <- phenodat[,c(1:8)]

phenodat$sex <- as.character(phenodat$sex)
phenodat$sex <- ifelse(phenodat$sex == '1', 'female', phenodat$sex)
phenodat$sex <- ifelse(phenodat$sex == '2', 'male', phenodat$sex)

phenodat_miss <- phenodat[!complete.cases(phenodat), ]
#12 missing! 

write_rds(phenodat, file='RNA_seq_pheno_data_checked.rds')
write_rds(phenodat_miss, file='RNA_seq_pheno_missing_data_checked.rds')

sampleinfo <- sampleinfo %>% select(SampleID, PAH)

#Check if the right data is ready for the design matrix
abc1 <- as.data.frame(colnames(counts))
colnames(abc1) <- 'SampleID'
abc2 <- as.data.frame(rownames(sampleinfo))
colnames(abc2) <- 'SampleID'
abc3 <- merge(abc1, abc2, by='SampleID')
#this returns with 453 rows, so perfect match. 

condition <- factor(sampleinfo$PAH)
levels(condition)

# NEW: Design matrix includes patient ID and condition
design <- model.matrix(~ condition)
rownames(design) <- rownames(sampleinfo)
view(design)

#======================#
# Filter lowly expressed genes
#======================#
counts <- counts %>% select(all_of(sampleinfo$SampleID))

myCPM <- cpm(counts)
#pick CPM of 0.5 and in at least 2 or more samples
thresh <- myCPM > 0.5
keep <- rowSums(thresh) >= 2
counts.keep <- counts[keep,]
dim(counts.keep)
#this retains 30,000 genes which may be a bit too relaxed? 

# Stricter filtering
thresh <- myCPM > 1
keep <- rowSums(thresh) >= 2
counts.keep <- counts[keep,]
dim(counts.keep)
#24600 genes kept, seems more reasonable 

# Use filtered counts for DGE
counts2 <- counts.keep

#======================#
# Differential Expression Analysis
#======================#
dge <- DGEList(counts = counts2)
dge <- calcNormFactors(dge)
#dge2 <- estimateDisp(dge, design)

#==================================================================
#add in PC & WBC adjustment
PCs_and_WBCs_RNA_seq$SampleID <- ifelse(PCs_and_WBCs_RNA_seq$SampleID == '1.0499999999999999E+81', 'OMITTED FOR PRIVACY', PCs_and_WBCs_RNA_seq$SampleID)

sampleinfo2 <- merge(sampleinfo, PCs_and_WBCs_RNA_seq, by = 'SampleID')
sampleinfo2_keep <- sampleinfo2 %>% group_by(External.ID) %>% filter(n() > 1) %>% slice(1) %>%  ungroup()
sampleinfo2_remove <- sampleinfo2 %>% group_by(External.ID) %>% filter(n() < 2) 
sampleinfo2 <- rbind(sampleinfo2_remove, sampleinfo2_keep)

sampleinfo2 <- sampleinfo2[,c(1:3,6:27)]


#update design matrix:
sampleinfo2 <- as.data.frame(sampleinfo2)
rownames(sampleinfo2) <- sampleinfo2$SampleID
sampleinfo2$PAH <- factor(sampleinfo2$PAH)
levels(sampleinfo2$PAH)
#first do for PCs Chris Rhodes adjusted too + Rowena's CIBER leukocytes
design <- model.matrix(~ PAH + ., data = sampleinfo2[, c(2,4:14,16,18,22:25)])
#rownames(design) <- rownames(sampleinfo2)
design_forleukos <- design

sample_for_clust <- sampleinfo2[,1:2]
write_rds(sample_for_clust, file='Clustering_samples_ready.rds')

#now do the DGE with adjustment
dge2 <- estimateDisp(dge, design)

#get log2CPM for actual expression
logCPM <- cpm(dge2, log = TRUE, prior.count = 1)

#Fit model and run test
fit <- glmQLFit(dge2, design)
test1 <- glmQLFTest(fit, coef='PAH1')
topTags(test1)

#Get 95% CI
logFC <- test1$table$logFC
Fstat <- test1$table$F
df <- test1$df.total
SE_logFC <- rep(NA, length(Fstat))
valid <- which(!is.na(Fstat) & Fstat > 0)
#7 variables dropped
SE_logFC[valid] <- abs(logFC[valid]) / sqrt(Fstat[valid])
names(SE_logFC) <- rownames(test1$table)

t_crit <- qt(0.975, df = df)

# Compute CIs
CI_lower <- logFC - t_crit * SE_logFC
CI_upper <- logFC + t_crit * SE_logFC

# Add to results
results_PAH <- as.data.frame(topTags(test1, n = 40000))
results_PAH$SE_logFC <- SE_logFC[rownames(results_PAH)]
results_PAH$logFC_lower <- results_PAH$logFC - t_crit * results_PAH$SE_logFC
results_PAH$logFC_upper <- results_PAH$logFC + t_crit * results_PAH$SE_logFC
write_rds(results_PAH, file='PAH_vs_control_PC_CIBER_adjusted_RNA_seq_11Now25.rds')
des_sum <- as.data.frame(design)
summary(as.factor(des_sum$PAH1))

IFN_genes <-c("OAS3", "ISG15", "RSAD2", "IFIT3", "IFI44", "OASL", "IFI44L", "IFIT1", "USP18", "XAF1")
BMP_genes <- c("PTGS2", "ARL4C", "ID2", "ID3", "NOTCH1", "NOTCH2", "SMAD1", "SMAD5")

results_PAH$gene <- rownames(results_PAH)
results_PAH$type <- ifelse(results_PAH$gene %in% IFN_genes, 'IFN', NA)
results_PAH$type <- ifelse(results_PAH$gene %in% BMP_genes, 'BMP', results_PAH$type)

res2 <- results_PAH %>% filter(!is.na(type))

gene_levels <- intersect(c(BMP_genes, IFN_genes), res2$gene)
res2$gene <- factor(res2$gene, levels = gene_levels)

ggplot(res2, aes(x = gene, y = logFC, colour = type)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin = logFC_lower, ymax = logFC_upper), width = 0.2) +
  xlab("Gene") +
  ylab("Log Fold Change (95% CI)") +
  ggtitle("Log Fold Change for BMP and IFN Genes") +
  theme_pubr() +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#only get IFN genes
ggplot(res2 %>% filter(type == 'IFN'), aes(x = gene, y = logFC, colour = type)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin = logFC_lower, ymax = logFC_upper), width = 0.2) +
  xlab("Gene") +
  ylab("Log Fold Change (95% CI)") +
  ggtitle("Log Fold Change for IFN Genes") +
  theme_pubr() +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

results_PAH <- results_PAH %>% filter(gene %in% genes_keep)
#214 kept (of 221 in dataset)
results_PAH$FDR2 <- p.adjust(results_PAH$PValue, method='fdr')

#plot volcano 
EnhancedVolcano(results_PAH,
                lab = rownames(results_PAH),
                x = 'logFC',
                y = 'FDR2',
                xlim = c(-1, 1),
                ylim=c(0,25),
                title = 'volcano plot of interferon transcript in cases vs controls',
                pCutoff = 0.05,
                FCcutoff = 0.1,
                pointSize = 3.0,
                labSize = 4.0)

EnhancedVolcano(results_PAH,
                lab = rownames(results_PAH),
                x = 'logFC',
                y = 'FDR2',
                xlim = c(-1, 1),
                ylim = c(0, 25),
                title = 'ISGs: PAH vs controls',
                subtitle = NULL,
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~ "FDR"),
                pCutoff = 0.05,
                FCcutoff = 0.1,
                pointSize = 3.0,
                labSize = 4.0,
                legendLabels = c('NS', expression(Log[2]~FC), 'FDR', expression(FDR~and~Log[2]~FC)),
                legendPosition = 'bottom')


#get heatmap for design matrix
cor_matrix <- cor(design, method = "spearman", use = "pairwise.complete.obs")

cor_matrix <- cor_matrix[-1,]
cor_matrix <- cor_matrix[,-1]


pheatmap(cor_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         display_numbers = FALSE,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         main = "Spearman Correlation (design matrix)",
         fontsize_row = 8,
         fontsize_col = 8)



res2$fdr2 <- p.adjust(res2$PValue, method='fdr')

###################################################################333
#now do the expression analysis
all_genes <- c(IFN_genes, BMP_genes)

genes_present <- intersect(all_genes, rownames(logCPM))
expr_data <- logCPM[genes_present, ]

sample_info <- data.frame(sample = colnames(expr_data),condition = design[, 'PAH1'])

expr_long <- expr_data %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "expression") %>%
  left_join(sample_info, by = "sample") %>%
  mutate(gene_type = case_when(gene %in% IFN_genes ~ "IFN",gene %in% BMP_genes ~ "BMP",
                               TRUE ~ "Other"),condition = factor(condition, levels = c(0, 1), labels = c("Control", "PAH")))


ggplot(expr_long, aes(x = gene, y = expression, fill = condition)) +
  geom_boxplot() +
  facet_wrap(~gene_type, scales = "free_x") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Gene",
    y = "Log2 CPM",
    title = "Expression of IFN and BMP Genes",
    fill = "Condition"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)


#==============================================================================
#heatmap ifn genes

#select the IFN genes of relevance
log2cpm <- as.data.frame(logCPM)
log2cpm <- log2cpm[rownames(log2cpm) %in% genes_keep, ]
log2cpm <- as.data.frame(t(log2cpm))
nameslog <- names(log2cpm)

#scale it
log2cpm$SampleID <- rownames(log2cpm)
sampl <- sampleinfo2 %>% select(SampleID, PAH)

log2cpm <- merge(log2cpm, sampl, by='SampleID')
write_rds(log2cpm, file='IFN_log_2_cpm.rds')

log2cpm[,2:215] <- lapply(log2cpm[,2:215], function(x) as.numeric(scale(x)))
log2cpm <- as.data.frame(log2cpm)

log2cpm$PAH <- as.character(log2cpm$PAH)
log2cpm$PAH <- ifelse(log2cpm$PAH == '1', 'PAH', 'Control')
log2cpm$PAH <- as.factor(log2cpm$PAH)

library(pheatmap)
annotation_df <- data.frame(PAH_Status = log2cpm$PAH)
rownames(annotation_df) <- log2cpm$SampleID

# Move SampleID back to row names so they appear on the heatmap
rownames(log2cpm) <- log2cpm$SampleID
heatmap_data <- log2cpm[, !(names(log2cpm) %in% c("SampleID", "PAH"))]
heatmap_matrix <- as.matrix(heatmap_data)
my_colors <- colorRampPalette(c("blue", "white", "red"))(100)
my_annotation_colors <- list(PAH_Status = c(PAH = "darkorange", Control = "purple") )
pheatmap(heatmap_matrix,
         color = my_colors,
         annotation_row = annotation_df,
         annotation_colors = my_annotation_colors, 
         cluster_rows = TRUE,   
         cluster_cols = TRUE,    
         show_rownames = FALSE,  
         show_colnames = FALSE,  
         fontsize_row = 8)      

#make complex heatmap for this for within group clustering
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)

row_ha <- rowAnnotation(
  Disease = log2cpm$PAH,
  col = list(PAH_Status = c(PAH = "darkorange", Control = "purple")))

Heatmap(heatmap_matrix, 
        name = "Z-scored log2CPM of IFN response genes", 
        col = my_colors,  
        right_annotation = row_ha,
        row_split = log2cpm$PAH, 
        cluster_rows = TRUE,     
        cluster_columns = TRUE,  
        show_row_names = FALSE,  
        show_column_names = FALSE, 
        row_gap = unit(2, "mm"),   
        row_title = NULL)          


disease_counts <- table(log2cpm$PAH)
split_labels <- paste0(log2cpm$PAH, " (n=", disease_counts[log2cpm$PAH], ")")


# 1. Create your heatmap keeping Cases and Controls separated
ht <- Heatmap(heatmap_matrix, 
              name = "Z-scored log2CPM of IFN response genes", 
              col = my_colors,  
              right_annotation = row_ha,
              row_split = log2cpm$PAH,   # <--- Keeps PAH and Control separated
              cluster_rows = TRUE,     
              cluster_columns = TRUE,  
              show_row_names = FALSE,  
              show_column_names = FALSE, 
              row_gap = unit(2, "mm"),   
              row_title_rot = 0) 

# 2. You MUST use draw() to calculate the clustering behind the scenes
ht_drawn <- draw(ht)

# 3. Extract the separate dendrograms
dends <- row_dend(ht_drawn)

# 4. Count how many people (leaves) are in the two main branches of each dendrogram
# (Using [[1]] for the first major branch and [[2]] for the second major branch)

# Control branches
control_branch_1 <- length(labels(dends[["Control"]][[1]]))
control_branch_2 <- length(labels(dends[["Control"]][[2]]))

# PAH branches
pah_branch_1 <- length(labels(dends[["PAH"]][[1]]))
pah_branch_2 <- length(labels(dends[["PAH"]][[2]]))

# 5. Print the counts!
cat("Control Group -> Branch 1:", control_branch_1, "| Branch 2:", control_branch_2, "\n")
cat("PAH Group     -> Branch 1:", pah_branch_1, "| Branch 2:", pah_branch_2, "\n")

#13/131 are high IFN in controls (9.9%)
#64/370 are high IFN in the cases (17.3%)

# Control branches
control_branch_1_ids <- labels(dends[["Control"]][[1]])
control_branch_2_ids <- labels(dends[["Control"]][[2]])

# PAH branches
pah_branch_1_ids <- labels(dends[["PAH"]][[1]])
pah_branch_2_ids <- labels(dends[["PAH"]][[2]])

# Print them
b1 <- cat("Control Branch 1 IDs:\n"); print(control_branch_1_ids)
cat("Control Branch 2 IDs:\n"); print(control_branch_2_ids)
cat("PAH Branch 1 IDs:\n");     print(pah_branch_1_ids)
cat("PAH Branch 2 IDs:\n");     print(pah_branch_2_ids)

# High IFN Controls (13)
high_ifn_controls <- c("OMITTED FOR PRIVACY")

# High IFN PAH (64)
high_ifn_pah <- c("OMITTED FOR PRIVACY")

high_ifn <- Complete_RNA_seq_cohort_pheno %>% filter(SampleID %in% high_ifn_controls | SampleID %in% high_ifn_pah)
high_ifn <- high_ifn %>% select(SampleID, CohortID, OpenClinicaID)
high_ifn <- merge(high_ifn, sampleinfo, by='SampleID')
getwd()
write_rds(high_ifn, file='Raised_IFN_groups_log2cpm.rds')
# ---------------------------------------------------------
# STEP 1: Calculate the dendrograms behind the scenes
# ---------------------------------------------------------
# We run your base heatmap just to let R do the math
ht_base <- Heatmap(heatmap_matrix, row_split = log2cpm$PAH)

# You must run draw() to generate the trees (a plot will appear, just ignore it for now)
ht_drawn <- draw(ht_base)
dend_list <- row_dend(ht_drawn)

# ---------------------------------------------------------
# STEP 2: Extract the branches and build custom labels
# ---------------------------------------------------------
# Create an empty vector to hold our new labels
new_splits <- character(nrow(heatmap_matrix))

# We need to map the patients safely, checking if your matrix has rownames
has_rownames <- !is.null(rownames(heatmap_matrix))
if(has_rownames) names(new_splits) <- rownames(heatmap_matrix)

# We will dynamically store the order so the final plot looks right
level_order <- c()

# Loop through your disease groups (Control, then PAH)
for(group in names(dend_list)) {
  
  # Get the exact row identities for the top and bottom branches
  branch_top <- labels(dend_list[[group]][[1]]) 
  branch_bot <- labels(dend_list[[group]][[2]]) 
  
  # Create the text you want printed on the plot
  lbl_top <- paste0(group, "\nTop Branch\n(n=", length(branch_top), ")")
  lbl_bot <- paste0(group, "\nBottom Branch\n(n=", length(branch_bot), ")")
  
  # Save the order for later
  level_order <- c(level_order, lbl_top, lbl_bot)
  
  # Assign the text to the correct rows
  if(has_rownames) {
    new_splits[branch_top] <- lbl_top
    new_splits[branch_bot] <- lbl_bot
  } else {
    new_splits[as.integer(branch_top)] <- lbl_top
    new_splits[as.integer(branch_bot)] <- lbl_bot
  }
}

# Turn it into a factor so the heatmap reads it in the correct top-to-bottom order
new_splits <- factor(new_splits, levels = level_order)

# ---------------------------------------------------------
# STEP 3: Draw the final visually labelled heatmap
# ---------------------------------------------------------
Heatmap(heatmap_matrix, 
        name = "Z-scored log2CPM of IFN response genes", 
        col = my_colors,  
        right_annotation = row_ha,
        row_split = new_splits,      # <--- We feed in our custom 4-block labels here!
        cluster_rows = TRUE,     
        cluster_columns = TRUE,  
        show_row_names = FALSE,  
        show_column_names = FALSE, 
        row_gap = unit(2, "mm"),
        row_title_rot = 0)


#########################################################################
#update plot:
# --- Step 1: Create IFN status variable ---
log2cpm$IFN_status <- ifelse(
  rownames(log2cpm) %in% c(high_ifn_controls, high_ifn_pah),
  "High",
  "Normal"
)
log2cpm$IFN_status <- factor(log2cpm$IFN_status, levels = c("Normal", "High"))

# --- Step 2: Update row_ha to include both annotations ---
row_ha <- rowAnnotation(
  PAH_Status = log2cpm$PAH,
  IFN_Status = log2cpm$IFN_status,
  col = list(
    PAH_Status = c(PAH = "darkorange", Control = "purple"),
    IFN_Status = c(High = "darkred", Normal = "darkblue")
  )
)

# --- Step 3: Build ht with updated row_ha ---
ht <- Heatmap(heatmap_matrix, 
              name = "Z-scored log2CPM of IFN response genes", 
              col = my_colors,  
              right_annotation = row_ha,
              row_split = log2cpm$PAH,   # <--- Keeps PAH and Control separated
              cluster_rows = TRUE,     
              cluster_columns = TRUE,  
              show_row_names = FALSE,  
              show_column_names = FALSE, 
              row_gap = unit(2, "mm"),   
              row_title_rot = 0) 

# --- Step 4: Draw as before ---
ht_drawn <- draw(ht)
dend_list <- row_dend(ht_drawn)

row_ha <- rowAnnotation(
  PAH_Status = log2cpm$PAH,
  IFN_Status = log2cpm$IFN_status,
  col = list(
    PAH_Status = c(PAH = "darkorange", Control = "purple"),
    IFN_Status = c(High = "darkred", Normal = "darkblue")
  ),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    PAH_Status = list(title = "PAH"),
    IFN_Status = list(title = "IFN")
  )
)

ht_final <- Heatmap(heatmap_matrix, 
                    name = "Z-scored log2CPM", 
                    col = my_colors,  
                    right_annotation = row_ha,
                    row_split = log2cpm$PAH,      # simple 2-group split: Control vs PAH
                    cluster_rows = TRUE,     
                    cluster_columns = TRUE,  
                    show_row_names = FALSE,  
                    show_column_names = FALSE, 
                    row_gap = unit(2, "mm"),
                    row_title_rot = 0)

draw(ht_final)

#get stats
result <- prop.test(x = c(11, 90), n = c(88, 365))
result

# --- Step 2: Capture it as a grob ---
heatmap_grob <- grid.grabExpr(draw(ht_final))

combined <- plot_grid(
  heatmap_grob, ifn_surv$plot,      # or km_plot$plot if using survminer
  ncol = 2,
  rel_widths = c(1.2, 1)
)
#================================================================
#get CV
#select the IFN genes of relevance
logunscale <- as.data.frame(logCPM)
logunscale <- logunscale[rownames(logunscale) %in% genes_keep, ]
logunscale <- as.data.frame(t(logunscale))
logunscale$SampleID <- rownames(logunscale)
sampl <- sampleinfo2 %>% select(SampleID, PAH)
logunscale <- merge(logunscale, sampl, by='SampleID')
logunscale$PAH <- as.character(logunscale$PAH)
logunscale$PAH <- ifelse(logunscale$PAH == '1', 'PAH', 'Control')
logunscale$PAH <- as.factor(logunscale$PAH)

cv_results <- logunscale %>% group_by(PAH) %>% summarise(across(where(is.numeric), 
                                                                ~ (sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE)) * 100)) %>%
  pivot_longer(cols = -PAH, 
               names_to = "Gene", 
               values_to = "CV") %>%
  pivot_wider(names_from = PAH, 
              values_from = CV, 
              names_prefix = "CV_")

library(ggpubr)
cv2 <- cv_results %>% pivot_longer(2:3, names_to = 'Group', values_to = 'CV')
ggplot(cv2, aes(x=Group, y=CV, fill=Group)) + 
  geom_boxplot() +
  stat_compare_means(method='wilcox.test') +
  theme_pubr() +
  xlab('Group') +
  ylab('Coefficient of variation') + 
  ggtitle('CV in patients and controls')

summary(cv_results$CV_Control)
summary(cv_results$CV_PAH)

cv_results$cv_difference <- cv_results$CV_PAH - cv_results$CV_Control

#compare leukocyte levels! 

leukos <- as.data.frame(design_forleukos)
leukos <- leukos %>% pivot_longer(7:19, names_to = 'Leukocyte', values_to = 'CIBER_value')
leukos$group <- ifelse(leukos$PAH1 == 1, 'PAH', 'control')

ggplot(leukos, aes(x=group, y=CIBER_value, fill=group)) +
  geom_boxplot() +
  #stat_compare_means(method='wilcox.test', comparisons = list(c("PAH", "control")), y.position = c(0.5, 0.5, 0.5, 0.5, 0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5))  +
  xlab('Group') +
  ylab('CIBER value') +
  theme_pubr() +
  facet_wrap(~Leukocyte, scales='free_y')

library(dplyr)

pvals_df <- leukos %>%
  group_by(Leukocyte) %>%
  summarise(p.value = wilcox.test(CIBER_value ~ group)$p.value) %>%
  mutate(label = ifelse(p.value < 0.001, "p < 0.001", paste0("p = ", signif(p.value, 2))))
pvals_df$fdr <- p.adjust(pvals_df$p.value, method='fdr')

#=========================================================================================
#Add in proteomics
proteins_interferone_genes <- readRDS("~/PhD/Projects/Interferone response/proteins_interferone_genes.rds")
IFN_ligands <- readRDS("C:/Users/Gebruiker/Downloads/IFN_ligands.rds")
Source_data_proteomics <- readRDS("Myfilepath/Omics/Source_data_proteomics.rds")
proteomics_easy_conversion_key <- readRDS("Myfilepath/Omics/proteomics_easy_conversion_key.rds")

#import data on the controls and cases
Cohort_control_info_proteomics <- read_excel("Myfilepath/Omics/Cohort control info_proteomics.xlsx")
follow_up_data <- readRDS("Myfilepath/Dataset for CYNAPSE/follow_up_data.rds")
clin_dat_at_diagnosis_clean <- readRDS("Myfilepath/Dataset for CYNAPSE/clin_dat_at_diagnosis_clean.rds")

#calculate the age and get sex for visit 1
age_sex <- clin_dat_at_diagnosis_clean %>% select(id, age_diagnosis, diagnosis_date, sex)
timediff <- follow_up_data %>% filter(visit == 1) %>% select(id, visit_date)
age_sex <- merge(age_sex, timediff, by='id', all=T)
age_sex$timediff <- as.numeric(age_sex$visit_date - age_sex$diagnosis_date)/365
summary(age_sex$timediff)
#negative must be 0
age_sex$timediff <- ifelse(age_sex$timediff <0, 0, age_sex$timediff)
age_sex$age_at_bleed <- age_sex$age_diagnosis + age_sex$timediff

age_sex <- age_sex %>% select(id, age_at_bleed, sex)
crtl <- Cohort_control_info_proteomics %>% select(OpenClinica.ID, age_diagnosis.checked, sex.checked)
names(crtl) <- c('id', 'age_at_bleed', 'sex')

demographics_proteomics <- rbind(age_sex, crtl)
write_rds(demographics_proteomics, file='demographics_of_proteomics.rds')

key <- proteomics_easy_conversion_key %>% filter(EntrezGeneSymbol %in% genes_keep)
response <- c(key$target_true)
ligands <- c(colnames(IFN_ligands[,3:15]), 'IFNA5')
#fix empty colnames
names(Source_data_proteomics)[is.na(names(Source_data_proteomics)) | names(Source_data_proteomics) == ""] <- paste0("unnamed_", seq_len(sum(is.na(names(Source_data_proteomics)) | names(Source_data_proteomics) == "")))
setdiff(response, names(Source_data_proteomics))
#one gene not recorded despite in conversion key --> fix by omitting
response <- response[response %in% names(Source_data_proteomics)]
prot <- Source_data_proteomics %>% filter(Event %in% c('PATIENT_1', 'CONTROL_1')) %>% select(Event, OpenClinica.ID, all_of(response), all_of(ligands))
prot <- unique(prot)
summary(as.factor(prot$Event))
#now calculate logFC and pval
prot[,3:95] <- lapply(prot[,3:95], as.numeric)
prot <- as.data.frame(prot)

#merge with demographics 
#571 patients, lets see how many have demographics
prot <- merge(prot, demographics_proteomics, by.x='OpenClinica.ID', by.y='id')
#566, so dropped 5 people without demographics

summary(prot$sex)
summary(prot$age_at_bleed)
#2 NAs for age at bleeding 

prot %>%
  pivot_longer(-c(Event, OpenClinica.ID, sex, age_at_bleed), names_to = "Protein", values_to = "Value") %>%
  group_by(Protein) %>%
  summarise(mean_control = mean(Value[Event == "CONTROL_1"], na.rm = TRUE)) %>%
  summarise(mean_of_means = mean(mean_control, na.rm = TRUE))
#so centre normalised to controls

#get difference
results_prot <- prot %>%
  pivot_longer(-c(Event, OpenClinica.ID, sex, age_at_bleed), names_to = "Protein", values_to = "Value") %>%
  group_by(Protein) %>%
  summarise(
    mean_patient = mean(Value[Event == "PATIENT_1"], na.rm = TRUE),
    mean_control = mean(Value[Event == "CONTROL_1"], na.rm = TRUE),
    Normalised_difference = mean_patient - mean_control,        # it's fine to call this log2FC
    p_value = t.test(Value ~ Event)$p.value,
    .groups = "drop")

results_prot$fdr <- p.adjust(results_prot$p_value, method='fdr')
results_prot <- as.data.frame(results_prot)
rownames(results_prot) <- results_prot$Protein

EnhancedVolcano(results_prot,
                lab = rownames(results_prot),
                x = 'Normalised_difference',
                y = 'fdr',
                xlim = c(-1, 2),
                ylim=c(0,24),
                title = 'volcano plot of interferon transcript in cases vs controls',
                pCutoff = 0.05,
                FCcutoff = 0.1,
                pointSize = 3.0,
                labSize = 4.0)

##########################
#get adjusted analysis:
library(limma)

# ---- 1. Set Event as a factor with CONTROL_1 as reference ----
prot_adj <- prot %>%
  filter(Event %in% c("CONTROL_1", "PATIENT_1")) %>%
  mutate(
    Event = factor(Event, levels = c("CONTROL_1", "PATIENT_1")),
    sex = factor(sex)
  )

levels(prot_adj$Event)   # should print "CONTROL_1" "PATIENT_1"
table(prot_adj$Event, prot_adj$sex)

# ---- 2. Complete cases for the covariates ----
# limma/model.matrix drops rows with NA covariates, which would misalign
# the matrix columns with the design, so remove them up front
prot_adj <- prot_adj %>%
  filter(!is.na(Event), !is.na(sex), !is.na(age_at_bleed))

# ---- 3. Expression matrix: proteins in rows, samples in columns ----
meta_cols <- c("Event", "OpenClinica.ID", "sex", "age_at_bleed")

mat <- prot_adj %>%
  select(-all_of(meta_cols)) %>%
  as.matrix() %>%
  t()

colnames(mat) <- prot_adj$OpenClinica.id   # optional, for readability

# ---- 4. Design matrix and limma fit ----
design <- model.matrix(~ Event + age_at_bleed + sex, data = prot_adj)
colnames(design)   # check the Event coefficient name

fit <- lmFit(mat, design)
fit <- eBayes(fit)

res_limma <- topTable(fit,
                      coef = "EventPATIENT_1",
                      number = Inf,
                      adjust.method = "fdr",
                      sort.by = "none")

# ---- 5. Volcano plot ----
# Data-driven axis limits
x_pad  <- 0.1
xlims  <- c(min(res_limma$logFC) - x_pad, max(res_limma$logFC) + x_pad)
ymax   <- max(-log10(res_limma$adj.P.Val), na.rm = TRUE)
ylims  <- c(0, ceiling(ymax * 1.1 * 2) / 2)   # ~10% headroom, rounded to nearest 0.5

EnhancedVolcano(res_limma,
                lab = rownames(res_limma),
                x = 'logFC',
                y = 'adj.P.Val',
                xlim = xlims,
                ylim = ylims,
                title = 'Cases vs controls (adjusted for age and sex)',
                subtitle = NULL,
                ylab = bquote(~-Log[10]~ "FDR"),
                legendLabels = c('NS', expression(Log[2]~FC), 'FDR', expression(FDR~and~Log[2]~FC)),
                pCutoff = 0.05,
                FCcutoff = 0.1,
                pointSize = 3.0,
                labSize = 4.0,
                titleLabSize = 14,
                drawConnectors = TRUE,
                widthConnectors = 0.4,
                max.overlaps = Inf,
                legendPosition = 'bottom')


###############################################################################################
#full DGE for IFN genes
#set wd
setwd('Myfilepath/Interferone response')

library(tidyverse)
library(survival)
library(survminer)
library(data.table)
library(purrr)
library(tibble)
library(DescTools)
library(boot)
library(gridExtra)
library(cowplot)
library(ggalluvial)
library(ggpubr)
library(naniar)

#load in the data
clean_comorbidity <- readRDS("Myfilepath/Dataset for CYNAPSE/clean_comorbidity.rds")
clean_medication <- readRDS("Myfilepath/Dataset for CYNAPSE/clean_medication.rds")
clin_dat_at_diagnosis_clean <- readRDS("Myfilepath/Dataset for CYNAPSE/clin_dat_at_diagnosis_clean.rds")
therapy_type_per_patient <- readRDS("~/PhD/Projects/Other (quick analyses)/Analyses for Joe/therapy_type_per_patient.rds")
hrct <- read.csv("~/PhD/Projects/Other (quick analyses)/Radiology data/code to be checked Justine/2. summarising hrct change/hrct.csv")
Cenus_1Dec_2023_UK_PAH_cohort <- readRDS("Myfilepath/Census December 2023/Cenus_1Dec_2023_UK_PAH_cohort.rds")

#construct mortality dataframe to merge everything with down the line
mort_df <-  clin_dat_at_diagnosis_clean
mort_df <- merge(mort_df, Cenus_1Dec_2023_UK_PAH_cohort, by='id')

mort_df$census_cause <- ifelse(mort_df$census_cause == 'DEATH', 'death', mort_df$census_cause)

#make sure sub_date is complete and format the data
mort_df$sub_date <- as.Date(mort_df$death_or_census_date)
mort_df$sub_date2 <- fifelse(is.na(mort_df$sub_date), ymd("2023-12-01"), mort_df$sub_date)
mort_df$surv_time <- (mort_df$sub_date2 - mort_df$diagnosis_date)
mort_df[c('surv_time', 'day')] <- str_split_fixed(mort_df$surv_time, ' ', 2)
mort_df$surv_time <- as.numeric(mort_df$surv_time)
mort_df$surv_time <- mort_df$surv_time/365.2422
mort_df$event <- ifelse(mort_df$census_cause == 'death', 1, 0)
mort_df$event <- ifelse(is.na(mort_df$event), 0, mort_df$event)
mort_df$diagnosis_verified <- as.factor(mort_df$diagnosis_verified)
mort_df$sex <- as.factor(mort_df$sex)
mort_df <- mort_df %>% filter(!is.na(surv_time))
#remove negative survival times (due to census date)
mort_df <- mort_df %>% filter(surv_time >=0)
mort_df <- unique(mort_df)

#censor at 15 years for cohort
mort_df$event <- ifelse(mort_df$surv_time >15 & mort_df$event == 1, 0, mort_df$event)
mort_df$surv_time <- ifelse(mort_df$surv_time >15, 15, mort_df$surv_time)

#load in data
Raised_IFN_groups_log2cpm <- readRDS("~/PhD/Projects/Interferone response/Raised_IFN_groups_log2cpm.rds")
Complete_RNA_seq_cohort_pheno <- readRDS("Myfilepath/Omics/Complete_RNA_seq_cohort_pheno.rds")
Samples_RNA_Seq_to_keep <- readRDS("Myfilepath/Interferone response/Samples_RNA_Seq_to_keep.rds")

has_rna <- c(unique(Samples_RNA_Seq_to_keep$OpenClinicaID))
mort_df <- mort_df %>% filter(id %in% has_rna)
#373 patients kept with RNA
mort_df$IFN <- ifelse(mort_df$id %in% Raised_IFN_groups_log2cpm$OpenClinicaID, 'raised', 'not raised')

mdf <- mort_df

sobj <- Surv(mdf$surv_time, mdf$event, type='right')
sfit <- survfit(sobj ~ IFN, data=mdf)

ifn_surv <- ggsurvplot(sfit, 
                       data = mdf, 
                       pval = TRUE,           
                       risk.table = TRUE,     
                       conf.int = TRUE, 
                       xlim=c(0,15),
                       break.x.by=2.5,
                       palette = c('darkblue', 'darkred'), 
                       legend.labs = c('Normal IFN', 'Raised IFN'),
                       title = "Survival for IFN level")
ifn_surv


#==============================================================================================================================
#develop a decision tree for the clusters
library(caret)
library(rpart)
library(rpart.plot)
library(glmnet)
library(pROC)
IFN_log_2_cpm <- readRDS("~/PhD/Projects/Interferone response/IFN_log_2_cpm.rds")

#==========================================================================
#Visualise cutpoint for Chris Rhodes

ifn_plot <- IFN_log_2_cpm %>% select(SampleID, PAH, IFIH1, WARS, HELB, PSMB8, TAP1, TRPM2)
ifn_plot <- ifn_plot %>% filter(PAH == 1)

# Yes/No flags per gene threshold
ifn_plot <- ifn_plot %>%
  mutate(
    IFIH1_prop  = ifelse(IFIH1 >= 6,   'yes', 'no'),
    WARS_prop   = ifelse(WARS  >= 7.7, 'yes', 'no'),
    HELB_prop   = ifelse(HELB  >= 6.3, 'yes', 'no'),
    PSMB8_prop  = ifelse(PSMB8 >= 7.7, 'yes', 'no'),
    TAP1_prop   = ifelse(TAP1  >= 8.4, 'yes', 'no'),
    TRPM2_prop  = ifelse(TRPM2 >= 2.4, 'yes', 'no')
  )

# Define the condition text for each gene
conditions <- data.frame(
  Gene      = c('IFIH1', 'WARS', 'HELB', 'PSMB8', 'TAP1', 'TRPM2'),
  Condition = c('IFIH1 >= 6', 'WARS >= 7.7', 'HELB >= 6.3',
                'PSMB8 >= 7.7', 'TAP1 >= 8.4', 'TRPM2 >= 2.4'),
  stringsAsFactors = FALSE
)

# Summary table: n and % "yes" only, per gene, with condition as a column
summary_table <- ifn_plot %>%
  select(ends_with('_prop')) %>%
  pivot_longer(everything(), names_to = 'Gene', values_to = 'Status') %>%
  mutate(Gene = sub('_prop', '', Gene)) %>%
  filter(Status == 'yes') %>%
  count(Gene, name = 'n') %>%
  mutate(Percent = round(100 * n / nrow(ifn_plot), 1)) %>%
  left_join(conditions, by = 'Gene') %>%
  select(Gene, Condition, n, Percent) %>%
  as.data.frame()

summary_table

#######################################################################

IFN_log_2_cpm$high_ifn <- ifelse(IFN_log_2_cpm$SampleID %in% Raised_IFN_groups_log2cpm$SampleID, 'High.IFN', 'Normal.IFN')
write_rds(IFN_log_2_cpm[,c(1,216,217)], 'IFN_groups_log2CPM_HC_based_June26.rds')


set.seed(42)

sample_df <- IFN_log_2_cpm

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: PREPARE DATA
# ══════════════════════════════════════════════════════════════════════════════

# Binary outcome as factor — adjust column name to match your data
model_data <- data.frame(
  cluster = factor(sample_df$high_ifn,
                   levels = c('High.IFN', 'Normal.IFN')),
  sample_df[, 2:215])

cat("Total samples:", nrow(model_data), "\n")
cat("Class distribution:\n")
print(table(model_data$cluster))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2: 75/25 STRATIFIED TRAIN/VALIDATION SPLIT
# (validation set is locked away — not touched until the very end)
# ══════════════════════════════════════════════════════════════════════════════

train_index <- createDataPartition(model_data$cluster,p= 0.75,list = FALSE)

train_set <- model_data[train_index, ]
val_set   <- model_data[-train_index, ]

cat("\nTraining samples:  ", nrow(train_set), "\n")  # ~375
cat("Validation samples:", nrow(val_set),   "\n")   # ~126
cat("\nClass balance in training set:\n")
print(table(train_set$cluster))
cat("Class balance in validation set:\n")
print(table(val_set$cluster))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3: LASSO ON TRAINING SET ONLY
# ══════════════════════════════════════════════════════════════════════════════

# Prepare matrices for glmnet (requires matrix input, not data frame)
x_train <- as.matrix(train_set[, -1])   # predictors only
y_train <- train_set$cluster            # outcome

# Cross-validated LASSO to find the best lambda (penalty strength)
# This CV is internal to LASSO only — separate from the tree's CV later
lasso_cv <- cv.glmnet(
  x         = x_train,
  y         = y_train,
  family    = "binomial",   # binary outcome
  alpha     = 1,            # 1 = LASSO, 0 = Ridge, 0.5 = ElasticNet
  nfolds    = 10,
  type.measure = "auc"      # optimise for AUC (better than deviance for imbalanced data)
)

plot(lasso_cv, main = "LASSO: CV AUC vs Lambda")

# Two standard lambda choices:
# lambda.min  = best CV AUC (can be aggressive, keeps more variables)
# lambda.1se  = simpler model within 1 SE of best (more conservative, fewer variables)
# For variable selection before a decision tree, lambda.1se is usually preferable
best_lambda <- lasso_cv$lambda.1se
cat("\nBest lambda (1se):", round(best_lambda, 5), "\n")

# Extract non-zero coefficients = selected transcripts
lasso_coefs    <- coef(lasso_cv, s = "lambda.1se")
selected_vars  <- rownames(lasso_coefs)[lasso_coefs[, 1] != 0]
selected_vars  <- selected_vars[selected_vars != "(Intercept)"]  # drop intercept

cat("Transcripts selected by LASSO:", length(selected_vars), "out of 214\n")
cat("\nSelected transcripts:\n")
print(selected_vars)

# Safety check — if LASSO selects 0 variables, relax to lambda.min
if (length(selected_vars) == 0) {
  cat("\nWARNING: lambda.1se selected 0 variables. Switching to lambda.min\n")
  best_lambda   <- lasso_cv$lambda.min
  lasso_coefs   <- coef(lasso_cv, s = "lambda.min")
  selected_vars <- rownames(lasso_coefs)[lasso_coefs[, 1] != 0]
  selected_vars <- selected_vars[selected_vars != "(Intercept)"]
  cat("Transcripts selected (lambda.min):", length(selected_vars), "\n")
  print(selected_vars)
}

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4: BUILD DECISION TREE ON LASSO-SELECTED VARIABLES WITH 10-FOLD CV
# ══════════════════════════════════════════════════════════════════════════════

# Subset training and validation sets to LASSO-selected variables only
train_lasso <- train_set[, c("cluster", selected_vars)]
val_lasso   <- val_set[,   c("cluster", selected_vars)]

# 10-fold CV control for the decision tree
tree_control <- trainControl(
  method          = "cv",
  number          = 10,
  savePredictions = "final",
  classProbs      = TRUE,
  summaryFunction = twoClassSummary   # reports ROC/AUC, Sens, Spec
)

# Train decision tree with CV to find best complexity parameter (CP)
m_tree <- train(
  cluster ~ .,
  data       = train_lasso,
  method     = "rpart",
  trControl  = tree_control,
  metric     = "ROC",          # optimise for AUC during CV
  tuneLength = 20              # tries 20 CP values and picks the best
)

cat("\n--- Decision Tree: 10-fold CV Results (on training set) ---\n")
print(m_tree)

# CV accuracy and AUC across folds
cat("\nBest CP value:", m_tree$bestTune$cp, "\n")
cat("CV AUC:       ", round(max(m_tree$results$ROC), 3), "\n")

# Plot CP tuning
plot(m_tree, main = "Decision Tree: CV AUC vs Complexity Parameter (CP)")

# Plot the final pruned tree
rpart.plot(
  m_tree$finalModel,
  type  = 4,
  extra = 104,   
  tweak = 1.2, # shows % of observations and class probability at each node
  main  = "Decision Tree — LASSO-selected transcripts"
)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5: VALIDATE ON HELD-OUT 25% (used for the first time here)
# ══════════════════════════════════════════════════════════════════════════════

val_preds <- predict(m_tree, newdata = val_lasso, type = "raw")
val_probs <- predict(m_tree, newdata = val_lasso, type = "prob")[, "High.IFN"]  # ← changed

cm <- confusionMatrix(
  data      = val_preds,
  reference = val_lasso$cluster,
  positive  = "High.IFN"   # ← changed
)

cat("\n========== VALIDATION SET RESULTS ==========\n")
print(cm)
cat("\nSummary:\n")
cat("Accuracy:   ", round(cm$overall["Accuracy"]    * 100, 1), "%\n")
cat("Sensitivity:", round(cm$byClass["Sensitivity"] * 100, 1), "%\n")
cat("Specificity:", round(cm$byClass["Specificity"] * 100, 1), "%\n")
cat("Kappa:      ", round(cm$overall["Kappa"],              3), "\n")

roc_val <- roc(
  response  = val_lasso$cluster,
  predictor = val_probs,
  levels    = c("Normal.IFN", "High.IFN")   # ← changed (negative first, positive second)
)

plot(
  roc_val,
  main = paste0("ROC Curve — Validation Set (AUC = ", round(auc(roc_val), 3), ")"),
  col  = "steelblue",
  lwd  = 2
)
abline(a = 0, b = 1, lty = 2, col = "grey")
cat("\nValidation AUC:", round(auc(roc_val), 3), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 6: SUMMARY OUTPUT
# ══════════════════════════════════════════════════════════════════════════════

cat("\n========== PIPELINE SUMMARY ==========\n")
cat("Total samples:          ", nrow(model_data),                          "\n")
cat("Training samples:       ", nrow(train_set),                           "\n")
cat("Validation samples:     ", nrow(val_set),                             "\n")
cat("Predictors (original):   214\n")
cat("Predictors (post-LASSO):", length(selected_vars),                     "\n")
cat("Decision tree CP:       ", round(m_tree$bestTune$cp, 5),              "\n")
cat("CV AUC (training):      ", round(max(m_tree$results$ROC), 3),         "\n")
cat("Validation AUC:         ", round(auc(roc_val), 3),                    "\n")
cat("Validation Accuracy:    ", round(cm$overall["Accuracy"] * 100, 1), "%\n")
cat("Selected transcripts:\n")
print(selected_vars)

#########################################################################################

#' ---
#' title: Visualisation of MOFA factors 
#' author: Eckart De Bie
#' output:
#'    html_document:
#'     toc: true
#'         
#' ---

setwd("Myfilepath/MOFA")
library(tidyverse)
library(ggpubr)

#load in the data
definitive_mofa_model_May25_PC_WBC_resids <- readRDS("~/PhD/Projects/MOFA/definitive_mofa_model_May25_PC_WBC_resids.rds")
definitive_mofa_model_factor_scores <- definitive_mofa_model_May25_PC_WBC_resids
MOFA_autoantibodies_weights_6Nov_25 <- read.csv2("~/PhD/Projects/MOFA/Supplement/MOFA_autoantibodies_weights_6Nov_25.csv")
MOFA_leukocytes_weights_6Nov_25 <- read.csv("~/PhD/Projects/MOFA/Supplement/MOFA_leukocytes_weights_6Nov_25.csv", sep=";")
MOFA_metabolomic_weights_6Nov_25 <- read.csv("~/PhD/Projects/MOFA/Supplement/MOFA_metabolomic_weights_6Nov_25.csv", sep=";")
MOFA_trascriptomic_weights_6Nov_25 <- read.csv2("~/PhD/Projects/MOFA/Supplement/MOFA_trascriptomic_weights_6Nov_25.csv")
MOFA_proteome_weights_6Nov_25 <- read.csv2("~/PhD/Projects/MOFA/Supplement/MOFA_proteome_weights_6Nov_25.csv")

#add the conversion here 
conversion_list_metabolites_3November26 <- readRDS("~/PhD/Projects/Interferone response/conversion_list_metabolites_3November26.rds")
MOFA_metabolomic_weights_6Nov_25 <- merge(MOFA_metabolomic_weights_6Nov_25, conversion_list_metabolites_3November26, by.x = 'X', by.y='namenew')
MOFA_metabolomic_weights_6Nov_25$X <- MOFA_metabolomic_weights_6Nov_25$Metabolite
MOFA_metabolomic_weights_6Nov_25 <- MOFA_metabolomic_weights_6Nov_25[,-c(17:20)]

#change rownames
# helper: replace values in df$X using a named vector (old = new)
recode_X <- function(df, map) {
  idx <- match(df$X, names(map))
  df$X[!is.na(idx)] <- unname(map[idx[!is.na(idx)]])
  df
}

# ---- Autoantibodies ----
# 1) generic suffix change: _norm -> _AAB
MOFA_autoantibodies_weights_6Nov_25$X <- sub("_norm$", "_AAB", MOFA_autoantibodies_weights_6Nov_25$X)

# 2) specific renames
MOFA_autoantibodies_weights_6Nov_25 <- recode_X(MOFA_autoantibodies_weights_6Nov_25, c(
  "Centromere.Protein.B..CENP.B._AAB"              = "Centromere Protein-B AAB",
  "Pyruvate.Dehydrogenase.from.porcine.hear_AAB"   = "Pyruvate Dehydrogenase AAB",
  "Scl.34..Fibrillarin.his.tag._AAB"               = "Scl-34 AAB",
  "Smith.antigen..Sm._AAB"                         = "Smith-antigen AAB",
  "Cardiolipin_AAB"              = "Cardiolipin AAB",
  "H2a.F2a2..and.H4.F2a1._AAB"   = "H2a (F2a2) and H4(F2a1) AAB",
  "Histone.type.IIA_AAB"               = "Histone type IIA AAB",
  "Jo.1_AAB"                         = "Jo-1 AAB",
  "La.SS.B.Antigen_AAB"              = "La/SS-B Antigen AAB",
  "Proteinase.3_AAB"   = "Proteinase-3 AAB",
  "Myeloperoxidase_AAB"               = "Myeloperoxidase AAB",
  "Mi.2b_AAB"                         = "Mi-2b AAB",
  "RNP.complex_AAB"              = "RNP-complex AAB",
  "Ro.SS.A..Antigen_AAB"   = "Ro/SS-A Antigen AAB",
  "SCL.70.Antigen_AAB"               = "Scl-70 Antigen AAB",
  "Thyroglobulin_AAB"                         = "Thyroglobulin AAB",
  "U1.snRNP.68_AAB"              = "U1-snRNP-68 AAB",
  "Thyroid.peroxidase_AAB"   = "Thyroid-peroxidase AAB",
  "Transglutaminase_AAB"               = "Transglutaminase AAB"
))

# ---- Leukocytes ----
MOFA_leukocytes_weights_6Nov_25 <- recode_X(MOFA_leukocytes_weights_6Nov_25, c(
  "ciber.B.cells.naive"                = "B.cells Naive",
  "ciber.B.cells.memory"               = "B-cells Memory",
  "ciber.Plasma.cells"                 = "Plasma cells",
  "ciber.T.cells.CD8"                  = "CD8 T-cells",
  "ciber.T.cells.CD4.naive"            = "CD4 T-cells Naive",
  "ciber.T.cells.CD4.memory.resting"   = "Memory CD4 T-cells resting",
  "ciber.T.cells.CD4.memory.activated" = "Memory CD4 T-cells activated",
  "ciber.T.cells.gamma.delta"          = "Gamma-delta T-cells",
  "ciber.NK.cells.resting"             = "NK cells resting",
  "ciber.NK.cells.activated"           = "NK cells activated",
  "ciber.Monocytes"                    = "Monocytes",
  "ciber.Macrophages.M0"               = "M0 Macrophages",
  "ciber.Macrophages.M1"               = "M1 Macrophages",
  "ciber.Dendritic.cells.resting"      = "Dendritic cells resting",
  "ciber.Dendritic.cells.activated"    = "Dendritic cells activated",
  "ciber.Mast.cells.resting"           = "Mast cells resting",
  "ciber.Eosinophils"                  = "Eosinophils",
  "ciber.Neutrophils"                  = "Neutrophils"
))

# check that nothing was missed
setdiff(MOFA_leukocytes_weights_6Nov_25$X, c("B.cells_Naive"))  
MOFA_leukocytes_weights_6Nov_25$X
MOFA_autoantibodies_weights_6Nov_25$X

#'Flag top factors for each modality and plot these
make_weight_plots <- function(df, factors, modality, n_top = 25) {
  
  df$X <- sub("_view_[0-9]+$", "", df$X)           
  n_use <- min(n_top, nrow(df))                    
  
  plots <- lapply(factors, function(f) {
    
    top_idx <- order(abs(df[[f]]), decreasing = TRUE)[seq_len(n_use)]
    
    plot_df <- df[top_idx, ] %>%
      transmute(Feature = X, Weight = .data[[f]]) %>%
      mutate(
        Feature = reorder(Feature, Weight),
        Direction = ifelse(Weight > 0, "Positive", "Negative"),
        label_x = ifelse(Weight > 0, -0.02, 0.02),
        label_h = ifelse(Weight > 0, 1, 0)
      )
    
    ggplot(plot_df, aes(x = Weight, y = Feature, colour = Direction)) +
      geom_vline(xintercept = 0, linetype = "dotted", colour = "grey30") +
      geom_segment(aes(x = 0, xend = Weight, yend = Feature), linewidth = 0.6) +
      geom_point(size = 2.5) +
      geom_text(aes(x = label_x, label = Feature, hjust = label_h),
                colour = "black", size = 3) +
      scale_colour_manual(values = c(Positive = "#B2182B", Negative = "#2166AC"),
                          guide = "none") +
      scale_x_continuous(breaks = seq(-1, 1, 0.5)) +
      coord_cartesian(xlim = c(-1, 1)) +
      labs(title = paste0(sub("Factor", "Factor ", f), " - top ", n_use, " ", modality, " weights"),
           x = "Weight", y = NULL) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major.y = element_blank(),
        axis.text.y        = element_blank(),
        axis.line.x        = element_line(colour = "black", linewidth = 0.5),
        legend.position    = "none"
      )
  })
  names(plots) <- factors
  plots
}

# ---- 2. Run for each modality ----
fac <- function(x) paste0("Factor", x)

modality_plots <- list(
  Transcriptomics = make_weight_plots(MOFA_trascriptomic_weights_6Nov_25,
                                      fac(c(1, 8, 9, 11, 12, 13, 14, 15)), "transcriptomic"),
  Proteomics      = make_weight_plots(MOFA_proteome_weights_6Nov_25,
                                      fac(c(2, 3, 4, 5, 6, 7, 10)), "proteomic"),
  Metabolomics    = make_weight_plots(MOFA_metabolomic_weights_6Nov_25,
                                      fac(c(2, 3, 5, 6, 10)), "metabolomic"),
  Autoantibodies  = make_weight_plots(MOFA_autoantibodies_weights_6Nov_25,
                                      fac(7), "autoantibody"),
  Leukocytes      = make_weight_plots(MOFA_leukocytes_weights_6Nov_25,
                                      fac(c(2, 3, 5, 6)), "leukocyte (CIBER)")
)

#'Make a figure
all_factors <- fac(1:15)

combined_list <- lapply(all_factors, function(f) {
  # collect this factor's plot from every modality that has it (others return NULL)
  p <- lapply(modality_plots, function(m) m[[f]])
  p <- Filter(Negate(is.null), p)
  if (length(p) == 0) return(NULL)
  ggarrange(plotlist = p, nrow = 1, ncol = length(p))
})
names(combined_list) <- all_factors
combined_list <- Filter(Negate(is.null), combined_list)

print(combined_list[["Factor2"]])
print(combined_list[["Factor3"]])
# Metabolomics, Factor 5
print(modality_plots$Metabolomics[["Factor5"]])

# Proteomics, Factor 5
print(modality_plots$Proteomics[["Factor5"]])
print(combined_list[["Factor10"]])
print(combined_list[["Factor6"]])

print(modality_plots$Metabolomics[["Factor6"]])

# Proteomics, Factor 5
print(modality_plots$Proteomics[["Factor6"]])

#' save these files 
for (f in names(combined_list)) {
  n_panels <- sum(sapply(modality_plots, function(m) !is.null(m[[f]])))
  ggsave(paste0("combined_weights_", f, ".png"), combined_list[[f]],
         width = 6 * n_panels, height = 7, dpi = 300, bg = "white")
}


#MOFA proteomics
setwd('Myfilepath/Interferone response')

#load libraries
library(tidyverse)
library(readxl)
library(data.table)
library(tidyr)
library(broom)
library(Publish)

definitive_mofa_model_factor_scores_adjused_for_filtering <- readRDS("~/PhD/Projects/MOFA/definitive_mofa_model_factor_scores_adjused_for_filtering.rds")

#get data for baseline table
definitive_mofa_model_factor_scores_adjused_for_filtering <- readRDS("~/PhD/Projects/MOFA/definitive_mofa_model_factor_scores_adjused_for_filtering.rds")
subsequent_REVEAL_Lite_worsening <- readRDS("Myfilepath/subsequent_REVEAL_Lite_worsening.rds")
incident_and_prevalent_cases <- readRDS("Myfilepath/incident_and_prevalent_cases.rds")
CCI_score_per_patient_complete <- readRDS("Myfilepath/CCI_score_per_patient_complete.rds")
clin_dat_at_diagnosis_clean <- readRDS("Myfilepath/Dataset for CYNAPSE/clin_dat_at_diagnosis_clean.rds")
therapy_type_per_patient <- readRDS("~/PhD/Projects/Other (quick analyses)/Analyses for Joe/therapy_type_per_patient.rds")
longitudional_score_differences <- readRDS("Myfilepath/longitudional_score_differences.rds")
follow_up_data <- readRDS("Myfilepath/Dataset for CYNAPSE/follow_up_data.rds")

clin <- clin_dat_at_diagnosis_clean %>% select(id, age_diagnosis, sex, diagnosis_verified, hb_pvr_calc, hb_pap_m, hb_rap_m, hb_pawp_m, cbt_card_ntprobnp_ngpl, fev1_fvc, lf_kco_pc, ep_1_type_6mwt, ep_1_distance_meters, functional_class, cbt_haem_hb_gpl, egfr_mdrd)
clin$mwd6 <- ifelse(clin$ep_1_type_6mwt == 'shuttle', clin$ep_1_distance_meters, NA)
rev1 <- longitudional_score_differences %>% filter(visit == 0) %>% select(id, reveal_lite)

clin <- merge(clin, rev1, by='id')
clin <- merge(clin, CCI_score_per_patient_complete, by='id')
clin <- merge(clin, incident_and_prevalent_cases, by='id')

clin <- clin %>% filter(id %in% definitive_mofa_model_factor_scores_adjused_for_filtering$samples)


time_df <- clin_dat_at_diagnosis_clean %>% select(id, diagnosis_date)
fudf <- follow_up_data %>% filter(visit == 1) %>% select(id, visit_date)
time_df <- merge(time_df, fudf, by='id')
time_df$visit_date <- as.Date(time_df$visit_date)
time_df$time_to_sampling <- as.numeric(time_df$visit_date - time_df$diagnosis_date)/365.24
time_df <- time_df %>% select(id, time_to_sampling)

clin <- merge(clin, time_df, by='id')
clin$diagnosis_verified <- droplevels(clin$diagnosis_verified)
clin <- merge(clin, ethnicity, by='id')
clin$diagnosis_verified <- as.factor(clin$diagnosis_verified)
clin$ethnic_category <- as.factor(clin$ethnic_category)
clin$functional_class <- as.factor(clin$functional_class)
clin$prevalent12months <- as.factor(clin$prevalent12months)


df <- summary(univariateTable(sex ~ diagnosis_verified + age_diagnosis + sex + ethnic_category + prevalent12months + CCI_adj_score + time_to_sampling + hb_pvr_calc + hb_pap_m + hb_rap_m + hb_pawp_m + cbt_card_ntprobnp_ngpl + fev1_fvc + lf_kco_pc + mwd6 + functional_class + reveal_lite + cbt_haem_hb_gpl + egfr_mdrd, data=clin, show.totals=TRUE))
write.csv2(df, file='MOFA_baseline.csv')


#add validation cohort
sample <- c(
  #PATIENT IDs OMITTED HERE FOR PRIVACY REASONS
)


clin <- clin_dat_at_diagnosis_clean %>% select(id, age_diagnosis, sex, diagnosis_verified, hb_pvr_calc, hb_pap_m, hb_rap_m, hb_pawp_m, cbt_card_ntprobnp_ngpl, fev1_fvc, lf_kco_pc, ep_1_type_6mwt, ep_1_distance_meters, functional_class, cbt_haem_hb_gpl, egfr_mdrd)
clin$mwd6 <- ifelse(clin$ep_1_type_6mwt == 'shuttle', clin$ep_1_distance_meters, NA)
rev1 <- longitudional_score_differences %>% filter(visit == 0) %>% select(id, reveal_lite)

clin <- merge(clin, rev1, by='id', all=T)
clin <- merge(clin, CCI_score_per_patient_complete, by='id', all=T)
clin <- merge(clin, incident_and_prevalent_cases, by='id', all=T)

clin <- clin %>% filter(id %in% sample)


time_df <- clin_dat_at_diagnosis_clean %>% select(id, diagnosis_date)
fudf <- follow_up_data %>% filter(visit == 1) %>% select(id, visit_date)
time_df <- merge(time_df, fudf, by='id', all=T)
time_df$visit_date <- as.Date(time_df$visit_date)
time_df$time_to_sampling <- as.numeric(time_df$visit_date - time_df$diagnosis_date)/365.24
time_df <- time_df %>% select(id, time_to_sampling)

clin <- merge(clin, time_df, by='id', all=T)
clin$diagnosis_verified <- droplevels(clin$diagnosis_verified)

clin <- merge(clin, ethnicity, by='id', all=T)
clin$diagnosis_verified <- as.factor(clin$diagnosis_verified)
clin$ethnic_category <- as.factor(clin$ethnic_category)
clin$functional_class <- as.factor(clin$functional_class)
clin$prevalent12months <- as.factor(clin$prevalent12months)

clin <- unique(clin)
clin <- clin %>% filter(!is.na(diagnosis_verified))

df <- summary(univariateTable(sex ~ diagnosis_verified + age_diagnosis + sex + ethnic_category + prevalent12months + CCI_adj_score + time_to_sampling + hb_pvr_calc + hb_pap_m + hb_rap_m + hb_pawp_m + cbt_card_ntprobnp_ngpl + fev1_fvc + lf_kco_pc + mwd6 + functional_class + reveal_lite + cbt_haem_hb_gpl + egfr_mdrd, data=clin, show.totals=TRUE))
write.csv2(df, file='MOFA_validation.csv')
