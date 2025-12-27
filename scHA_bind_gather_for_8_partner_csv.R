library(wlab.block)
library(dplyr)
wt_aa<-"TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
#K55
## caculate interchain distance
interchain_distance_K55 <- bio3D_minimum_interchain_distances(
  input_file = "5mla.pdb",
  chain_query = "A",
  chain_target = "B")
colnames(interchain_distance_K55)[2:3]<-paste0(colnames(interchain_distance_K55)[2:3],"_K55")

## caculate distance from chain to a specific HETATM
chain_hetatm_dis_K55 <- bio3D_chain_hetatm_distances(
  input_file = "5mla.pdb",
  chain_query = "A",
  hetatm_target = "GNP"
)
colnames(chain_hetatm_dis_K55)[2:3]<-paste0(colnames(chain_hetatm_dis_K55)[2:3],"_GNP")
chain_hetatm_dis_K55 <- chain_hetatm_dis_K55 %>%
  rename(
    GXPMG_HAmin_ligand = HAmin_ligand_GNP,
    GXPMG_scHAmin_ligand = scHAmin_ligand_GNP
  )

colnames(chain_hetatm_dis_K55)[2:3]<-paste0(colnames(chain_hetatm_dis_K55)[2:3],"_K55")


anno_K55<-merge(interchain_distance_K55,chain_hetatm_dis_K55,by='Pos')
anno_K55[, K55_bind := 'no']
anno_K55[, K55_GNP_bind := 'no']
anno_K55[scHAmin_ligand_K55 <5 , K55_bind := 'yes']
anno_K55[GXPMG_scHAmin_ligand_K55 <5 , K55_GNP_bind := 'yes']



#####################
#####K13
## caculate interchain distance
interchain_distance_K13 <- bio3D_minimum_interchain_distances(
  input_file = "6h46.pdb",
  chain_query = "A",
  chain_target = "B")
colnames(interchain_distance_K13)[2:3]<-paste0(colnames(interchain_distance_K13)[2:3],"_K13")

## caculate distance from chain to a specific HETATM
chain_hetatm_dis_K13 <- bio3D_chain_hetatm_distances(
  input_file = "6h46.pdb",
  chain_query = "A",
  hetatm_target = "GNP"
)
colnames(chain_hetatm_dis_K13)[2:3]<-paste0(colnames(chain_hetatm_dis_K13)[2:3],"_GNP")
chain_hetatm_dis_K13 <- chain_hetatm_dis_K13 %>%
  rename(
    GXPMG_HAmin_ligand = HAmin_ligand_GNP,
    GXPMG_scHAmin_ligand = scHAmin_ligand_GNP
  )

colnames(chain_hetatm_dis_K13)[2:3]<-paste0(colnames(chain_hetatm_dis_K13)[2:3],"_K13")


anno_K13<-merge(interchain_distance_K13,chain_hetatm_dis_K13,by='Pos')
anno_K13[, K13_bind := 'no']
anno_K13[, K13_GNP_bind := 'no']
anno_K13[scHAmin_ligand_K13 <5 , K13_bind := 'yes']
anno_K13[GXPMG_scHAmin_ligand_K13 <5 , K13_GNP_bind := 'yes']



#############
###K19
## caculate interchain distance
interchain_distance_K19 <- bio3D_minimum_interchain_distances(
  input_file = "6h47.pdb",
  chain_query = "A",
  chain_target = "B")
colnames(interchain_distance_K19)[2:3]<-paste0(colnames(interchain_distance_K19)[2:3],"_K19")

## caculate distance from chain to a specific HETATM
chain_hetatm_dis_K19 <- bio3D_chain_hetatm_distances(
  input_file = "6h46.pdb",
  chain_query = "A",
  hetatm_target = "GNP"
)
colnames(chain_hetatm_dis_K19)[2:3]<-paste0(colnames(chain_hetatm_dis_K19)[2:3],"_GNP")
chain_hetatm_dis_K19 <- chain_hetatm_dis_K19 %>%
  rename(
    GXPMG_HAmin_ligand = HAmin_ligand_GNP,
    GXPMG_scHAmin_ligand = scHAmin_ligand_GNP
  )

colnames(chain_hetatm_dis_K19)[2:3]<-paste0(colnames(chain_hetatm_dis_K19)[2:3],"_K19")


anno_K19<-merge(interchain_distance_K19,chain_hetatm_dis_K19,by='Pos')
anno_K19[, K19_bind := 'no']
anno_K19[, K19_GNP_bind := 'no']
anno_K19[scHAmin_ligand_K19 <5 , K19_bind := 'yes']
anno_K19[GXPMG_scHAmin_ligand_K19 <5 , K19_GNP_bind := 'yes']



############################
##### RAF1
## caculate interchain distance
interchain_distance_RAF1 <- bio3D_minimum_interchain_distances(
  input_file = "6vjj.pdb",
  chain_query = "A",
  chain_target = "B")
colnames(interchain_distance_RAF1)[2:3]<-paste0(colnames(interchain_distance_RAF1)[2:3],"_RAF1")

## caculate distance from chain to a specific HETATM
chain_hetatm_dis_RAF1 <- bio3D_chain_hetatm_distances(
  input_file = "6vjj.pdb",
  chain_query = "A",
  hetatm_target = "GNP"
)
colnames(chain_hetatm_dis_RAF1)[2:3]<-paste0(colnames(chain_hetatm_dis_RAF1)[2:3],"_GNP")
chain_hetatm_dis_RAF1 <- chain_hetatm_dis_RAF1 %>%
  rename(
    GXPMG_HAmin_ligand = HAmin_ligand_GNP,
    GXPMG_scHAmin_ligand = scHAmin_ligand_GNP
  )

colnames(chain_hetatm_dis_RAF1)[2:3]<-paste0(colnames(chain_hetatm_dis_RAF1)[2:3],"_RAF1")


anno_RAF1<-merge(interchain_distance_RAF1,chain_hetatm_dis_RAF1,by='Pos')
anno_RAF1[, RAF1_bind := 'no']
anno_RAF1[, RAF1_GNP_bind := 'no']
anno_RAF1[scHAmin_ligand_RAF1 <5 , RAF1_bind := 'yes']
anno_RAF1[GXPMG_scHAmin_ligand_RAF1 <5 , RAF1_GNP_bind := 'yes']



###################
##K27
## caculate interchain distance
interchain_distance_K27 <- bio3D_minimum_interchain_distances(
  input_file = "5mlb.pdb",
  chain_query = "A",
  chain_target = "B")
colnames(interchain_distance_K27)[2:3]<-paste0(colnames(interchain_distance_K27)[2:3],"_K27")

## caculate distance from chain to a specific HETATM
chain_hetatm_dis_K27 <- bio3D_chain_hetatm_distances(
  input_file = "5mla.pdb",
  chain_query = "A",
  hetatm_target = "GNP"
)
colnames(chain_hetatm_dis_K27)[2:3]<-paste0(colnames(chain_hetatm_dis_K27)[2:3],"_GNP")
chain_hetatm_dis_K27 <- chain_hetatm_dis_K27 %>%
  rename(
    GXPMG_HAmin_ligand = HAmin_ligand_GNP,
    GXPMG_scHAmin_ligand = scHAmin_ligand_GNP
  )

colnames(chain_hetatm_dis_K27)[2:3]<-paste0(colnames(chain_hetatm_dis_K27)[2:3],"_K27")


anno_K27<-merge(interchain_distance_K27,chain_hetatm_dis_K27,by='Pos')
anno_K27[, K27_bind := 'no']
anno_K27[, K27_GNP_bind := 'no']
anno_K27[scHAmin_ligand_K27 <5 , K27_bind := 'yes']
anno_K27[GXPMG_scHAmin_ligand_K27 <5 , K27_GNP_bind := 'yes']




###################
##PI3KCG
## caculate interchain distance
interchain_distance_PI3KCG <- bio3D_minimum_interchain_distances(
  input_file = "1he8.pdb",
  chain_query = "B",
  chain_target = "A")
colnames(interchain_distance_PI3KCG)[2:3]<-paste0(colnames(interchain_distance_PI3KCG)[2:3],"_PI3KCG")

## caculate distance from chain to a specific HETATM
chain_hetatm_dis_PI3KCG <- bio3D_chain_hetatm_distances(
  input_file = "6vjj.pdb",
  chain_query = "A",
  hetatm_target = "GNP"
)
colnames(chain_hetatm_dis_PI3KCG)[2:3]<-paste0(colnames(chain_hetatm_dis_PI3KCG)[2:3],"_GNP")
chain_hetatm_dis_PI3KCG <- chain_hetatm_dis_PI3KCG %>%
  rename(
    GXPMG_HAmin_ligand = HAmin_ligand_GNP,
    GXPMG_scHAmin_ligand = scHAmin_ligand_GNP
  )

colnames(chain_hetatm_dis_PI3KCG)[2:3]<-paste0(colnames(chain_hetatm_dis_PI3KCG)[2:3],"_PI3KCG")

anno_PI3KCG<-merge(interchain_distance_PI3KCG,chain_hetatm_dis_PI3KCG,by='Pos')
anno_PI3KCG[, PI3KCG_bind := 'no']
anno_PI3KCG[, PI3KCG_GNP_bind := 'no']
anno_PI3KCG[scHAmin_ligand_PI3KCG <5 , PI3KCG_bind := 'yes']
anno_PI3KCG[GXPMG_scHAmin_ligand_PI3KCG <5 , PI3KCG_GNP_bind := 'yes']


###################这个计算不出来
##RALGDS
## caculate interchain distance
interchain_distance_RALGDS <- bio3D_minimum_interchain_distances(
  input_file = "1lfd.pdb",
  chain_query = "C",
  chain_target = "D")
colnames(interchain_distance_RALGDS)[2:3]<-paste0(colnames(interchain_distance_RALGDS)[2:3],"_RALGDS")

## caculate distance from chain to a specific HETATM
chain_hetatm_dis_RALGDS <- bio3D_chain_hetatm_distances(
  input_file = "6vjj.pdb",
  chain_query = "A",
  hetatm_target = "GNP"
)
colnames(chain_hetatm_dis_RALGDS)[2:3]<-paste0(colnames(chain_hetatm_dis_RALGDS)[2:3],"_GNP")
chain_hetatm_dis_RALGDS <- chain_hetatm_dis_RALGDS %>%
  rename(
    GXPMG_HAmin_ligand = HAmin_ligand_GNP,
    GXPMG_scHAmin_ligand = scHAmin_ligand_GNP
  )

colnames(chain_hetatm_dis_RALGDS)[2:3]<-paste0(colnames(chain_hetatm_dis_RALGDS)[2:3],"_RALGDS")

anno_RALGDS<-merge(interchain_distance_RALGDS,chain_hetatm_dis_RALGDS,by='Pos')
anno_RALGDS[, RALGDS_bind := 'no']
anno_RALGDS[, RALGDS_GNP_bind := 'no']
anno_RALGDS[scHAmin_ligand_RALGDS <5 , RALGDS_bind := 'yes']
anno_RALGDS[GXPMG_scHAmin_ligand_RALGDS <5 , RALGDS_GNP_bind := 'yes']




###################
##SOS1
## caculate interchain distance
interchain_distance_SOS1 <- bio3D_minimum_interchain_distances(
  input_file = "1nvw.pdb",
  chain_query = "Q",
  chain_target = "S")
colnames(interchain_distance_SOS1)[2:3]<-paste0(colnames(interchain_distance_SOS1)[2:3],"_SOS1")

## caculate distance from chain to a specific HETATM
chain_hetatm_dis_SOS1 <- bio3D_chain_hetatm_distances(
  input_file = "1nvw.pdb",
  chain_query = "Q",
  hetatm_target = "GNP"
)
colnames(chain_hetatm_dis_SOS1)[2:3]<-paste0(colnames(chain_hetatm_dis_SOS1)[2:3],"_GNP")
chain_hetatm_dis_SOS1 <- chain_hetatm_dis_SOS1 %>%
  rename(
    GXPMG_HAmin_ligand = HAmin_ligand_GNP,
    GXPMG_scHAmin_ligand = scHAmin_ligand_GNP
  )

colnames(chain_hetatm_dis_SOS1)[2:3]<-paste0(colnames(chain_hetatm_dis_SOS1)[2:3],"_SOS1")

anno_SOS1<-merge(interchain_distance_SOS1,chain_hetatm_dis_SOS1,by='Pos')
anno_SOS1[, SOS1_bind := 'no']
anno_SOS1[, SOS1_GNP_bind := 'no']
anno_SOS1[scHAmin_ligand_SOS1 <5 , SOS1_bind := 'yes']
anno_SOS1[GXPMG_scHAmin_ligand_SOS1 <5 , SOS1_GNP_bind := 'yes']





anno_final<-merge(anno_K13,anno_K19,by="Pos",all = T)
anno_final<-merge(anno_final,anno_K27,by="Pos",all = T)
anno_final<-merge(anno_final,anno_K55,by="Pos",all = T)
anno_final<-merge(anno_final,anno_RAF1,by="Pos",all = T)
anno_final<-merge(anno_final,anno_SOS1,by="Pos",all = T)
anno_final<-merge(anno_final,anno_RALGDS,by="Pos",all = T)
anno_final<-merge(anno_final,anno_PI3KCG,by="Pos",all = T)
write.csv(anno_final,"anno_final_for_8.csv")
