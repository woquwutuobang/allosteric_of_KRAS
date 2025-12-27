library(wlab.block)
library(dplyr)
wt_aa<-"TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"




######test 
interchain_distance <- bio3D_minimum_interchain_distances(
  input_file = "./7l0g.pdb",
  chain_query = "A",
  chain_target = "C")
colnames(interchain_distance)[2:3]<-paste0(colnames(interchain_distance)[2:3],"_VH")

colnames(interchain_distance)[1] <- c("Pos_real")
anno<-interchain_distance
anno[, VH_bind := 'no']
anno[scHAmin_ligand_VH <5 , VH_bind := 'yes']
anno[VH_bind=="yes",site_type:="Binding interface site"]
names(anno)
library(dplyr)
binding_positions <- anno %>%
  filter(site_type == "Binding interface site") %>%
  pull(Pos_real)
print(binding_positions)




#K55
## caculate interchain distance
interchain_distance_K55 <- bio3D_minimum_interchain_distances(
  input_file = "5mla.pdb",
  chain_query = "A",
  chain_target = "B")
colnames(interchain_distance_K55)[2:3]<-paste0(colnames(interchain_distance_K55)[2:3],"_K55")

## caculate distance from chain to a specific HETATM
chain_hetatm_dis_K55 <- bio3D_chain_hetatm_distances(
  input_file = "6vjj.pdb",
  chain_query = "A",
  hetatm_target = "GNP"
)
colnames(chain_hetatm_dis_K55)[2:3]<-paste0(colnames(chain_hetatm_dis_K55)[2:3],"_GTP")
chain_hetatm_dis_K55 <- chain_hetatm_dis_K55 %>%
  rename(
    GXPMG_HAmin_ligand = HAmin_ligand_GTP,
    GXPMG_scHAmin_ligand = scHAmin_ligand_GTP
  )

colnames(chain_hetatm_dis_K55)[2:3]<-paste0(colnames(chain_hetatm_dis_K55)[2:3],"_K55")


anno_K55<-merge(interchain_distance_K55,chain_hetatm_dis_K55,by='Pos')
anno_K55[, K55_bind := 'no']
anno_K55[, K55_GTP_bind := 'no']
anno_K55[scHAmin_ligand_K55 <5 , K55_bind := 'yes']
anno_K55[GXPMG_scHAmin_ligand_K55 <5 , K55_GTP_bind := 'yes']



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
  input_file = "6vjj.pdb",
  chain_query = "A",
  hetatm_target = "GNP"
)
colnames(chain_hetatm_dis_K13)[2:3]<-paste0(colnames(chain_hetatm_dis_K13)[2:3],"_GTP")
chain_hetatm_dis_K13 <- chain_hetatm_dis_K13 %>%
  rename(
    GXPMG_HAmin_ligand = HAmin_ligand_GTP,
    GXPMG_scHAmin_ligand = scHAmin_ligand_GTP
  )

colnames(chain_hetatm_dis_K13)[2:3]<-paste0(colnames(chain_hetatm_dis_K13)[2:3],"_K13")


anno_K13<-merge(interchain_distance_K13,chain_hetatm_dis_K13,by='Pos')
anno_K13[, K13_bind := 'no']
anno_K13[, K13_GTP_bind := 'no']
anno_K13[scHAmin_ligand_K13 <5 , K13_bind := 'yes']
anno_K13[GXPMG_scHAmin_ligand_K13 <5 , K13_GTP_bind := 'yes']



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
  input_file = "6vjj.pdb",
  chain_query = "A",
  hetatm_target = "GNP"
)
colnames(chain_hetatm_dis_K19)[2:3]<-paste0(colnames(chain_hetatm_dis_K19)[2:3],"_GTP")
chain_hetatm_dis_K19 <- chain_hetatm_dis_K19 %>%
  rename(
    GXPMG_HAmin_ligand = HAmin_ligand_GTP,
    GXPMG_scHAmin_ligand = scHAmin_ligand_GTP
  )

colnames(chain_hetatm_dis_K19)[2:3]<-paste0(colnames(chain_hetatm_dis_K19)[2:3],"_K19")


anno_K19<-merge(interchain_distance_K19,chain_hetatm_dis_K19,by='Pos')
anno_K19[, K19_bind := 'no']
anno_K19[, K19_GTP_bind := 'no']
anno_K19[scHAmin_ligand_K19 <5 , K19_bind := 'yes']
anno_K19[GXPMG_scHAmin_ligand_K19 <5 , K19_GTP_bind := 'yes']



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
colnames(chain_hetatm_dis_RAF1)[2:3]<-paste0(colnames(chain_hetatm_dis_RAF1)[2:3],"_GTP")
chain_hetatm_dis_RAF1 <- chain_hetatm_dis_RAF1 %>%
  rename(
    GXPMG_HAmin_ligand = HAmin_ligand_GTP,
    GXPMG_scHAmin_ligand = scHAmin_ligand_GTP
  )

colnames(chain_hetatm_dis_RAF1)[2:3]<-paste0(colnames(chain_hetatm_dis_RAF1)[2:3],"_RAF1")


anno_RAF1<-merge(interchain_distance_RAF1,chain_hetatm_dis_RAF1,by='Pos')
anno_RAF1[, RAF1_bind := 'no']
anno_RAF1[, RAF1_GTP_bind := 'no']
anno_RAF1[scHAmin_ligand_RAF1 <5 , RAF1_bind := 'yes']
anno_RAF1[GXPMG_scHAmin_ligand_RAF1 <5 , RAF1_GTP_bind := 'yes']



###################
##K27
## caculate interchain distance
interchain_distance_K27 <- bio3D_minimum_interchain_distances(
  input_file = "5o2s.pdb",
  chain_query = "A",
  chain_target = "B")
colnames(interchain_distance_K27)[2:3]<-paste0(colnames(interchain_distance_K27)[2:3],"_K27")

## caculate distance from chain to a specific HETATM
chain_hetatm_dis_K27 <- bio3D_chain_hetatm_distances(
  input_file = "6vjj.pdb",
  chain_query = "A",
  hetatm_target = "GNP"
)
colnames(chain_hetatm_dis_K27)[2:3]<-paste0(colnames(chain_hetatm_dis_K27)[2:3],"_GTP")
chain_hetatm_dis_K27 <- chain_hetatm_dis_K27 %>%
  rename(
    GXPMG_HAmin_ligand = HAmin_ligand_GTP,
    GXPMG_scHAmin_ligand = scHAmin_ligand_GTP
  )

colnames(chain_hetatm_dis_K27)[2:3]<-paste0(colnames(chain_hetatm_dis_K27)[2:3],"_K27")


anno_K27<-merge(interchain_distance_K27,chain_hetatm_dis_K27,by='Pos')
anno_K27[, K27_bind := 'no']
anno_K27[, K27_GTP_bind := 'no']
anno_K27[scHAmin_ligand_K27 <5 , K27_bind := 'yes']
anno_K27[GXPMG_scHAmin_ligand_K27 <5 , K27_GTP_bind := 'yes']


anno_final<-merge(anno_K13,anno_K19,by="Pos",all = T)
anno_final<-merge(anno_final,anno_K27,by="Pos",all = T)
anno_final<-merge(anno_final,anno_K55,by="Pos",all = T)
anno_final<-merge(anno_final,anno_RAF1,by="Pos",all = T)
write.csv(anno_final,"anno_final_ma.csv")

#############################
########dimer1 ##这个是α4-α5/α4-α5 dimer以A到B的距离
## caculate interchain distance
interchain_distance_dimer1 <- bio3D_minimum_interchain_distances(
  input_file = "7rsc.pdb",
  chain_query = "A",
  chain_target = "B")
colnames(interchain_distance_dimer1)[2:3]<-paste0(colnames(interchain_distance_dimer1)[2:3],"_dimer1")

## caculate distance from chain to a specific HETATM
chain_hetatm_dis_dimer1 <- bio3D_chain_hetatm_distances(
  input_file = "6vjj.pdb",
  chain_query = "A",
  hetatm_target = "GNP"
)
colnames(chain_hetatm_dis_dimer1)[2:3]<-paste0(colnames(chain_hetatm_dis_dimer1)[2:3],"_GTP")
chain_hetatm_dis_dimer1 <- chain_hetatm_dis_dimer1 %>%
  rename(
    GXPMG_HAmin_ligand = HAmin_ligand_GTP,
    GXPMG_scHAmin_ligand = scHAmin_ligand_GTP
  )

colnames(chain_hetatm_dis_dimer1)[2:3]<-paste0(colnames(chain_hetatm_dis_dimer1)[2:3],"_dimer1")


anno_dimer1<-merge(interchain_distance_dimer1,chain_hetatm_dis_dimer1,by='Pos')
anno_dimer1[, dimer1_bind := 'no']
anno_dimer1[, dimer1_GTP_bind := 'no']
anno_dimer1[scHAmin_ligand_dimer1 <5 , dimer1_bind := 'yes']
anno_dimer1[GXPMG_scHAmin_ligand_dimer1 <5 , dimer1_GTP_bind := 'yes']




#############################
########dimer2 ##这个是α4-α5/α4-α5 dimer以B到A的距离
## caculate interchain distance
interchain_distance_dimer2 <- bio3D_minimum_interchain_distances(
  input_file = "7rsc.pdb",
  chain_query = "B",
  chain_target = "A")
colnames(interchain_distance_dimer2)[2:3]<-paste0(colnames(interchain_distance_dimer2)[2:3],"_dimer2")

## caculate distance from chain to a specific HETATM
chain_hetatm_dis_dimer2 <- bio3D_chain_hetatm_distances(
  input_file = "6vjj.pdb",
  chain_query = "A",
  hetatm_target = "GNP"
)
colnames(chain_hetatm_dis_dimer2)[2:3]<-paste0(colnames(chain_hetatm_dis_dimer2)[2:3],"_GTP")
chain_hetatm_dis_dimer2 <- chain_hetatm_dis_dimer2 %>%
  rename(
    GXPMG_HAmin_ligand = HAmin_ligand_GTP,
    GXPMG_scHAmin_ligand = scHAmin_ligand_GTP
  )

colnames(chain_hetatm_dis_dimer2)[2:3]<-paste0(colnames(chain_hetatm_dis_dimer2)[2:3],"_dimer2")


anno_dimer2<-merge(interchain_distance_dimer2,chain_hetatm_dis_dimer2,by='Pos')
anno_dimer2[, dimer2_bind := 'no']
anno_dimer2[, dimer2_GTP_bind := 'no']
anno_dimer2[scHAmin_ligand_dimer2 <5 , dimer2_bind := 'yes']
anno_dimer2[GXPMG_scHAmin_ligand_dimer2 <5 , dimer2_GTP_bind := 'yes']






#############################
########dimer4 ##这个是α4-α5/β dimer以B到A的距离
## caculate interchain distance
interchain_distance_dimer4 <- bio3D_minimum_interchain_distances(
  input_file = "7rse.pdb",
  chain_query = "B",
  chain_target = "A")
colnames(interchain_distance_dimer4)[2:3]<-paste0(colnames(interchain_distance_dimer4)[2:3],"_dimer4")

## caculate distance from chain to a specific HETATM
chain_hetatm_dis_dimer4 <- bio3D_chain_hetatm_distances(
  input_file = "6vjj.pdb",
  chain_query = "A",
  hetatm_target = "GNP"
)
colnames(chain_hetatm_dis_dimer4)[2:3]<-paste0(colnames(chain_hetatm_dis_dimer4)[2:3],"_GTP")
chain_hetatm_dis_dimer4 <- chain_hetatm_dis_dimer4 %>%
  rename(
    GXPMG_HAmin_ligand = HAmin_ligand_GTP,
    GXPMG_scHAmin_ligand = scHAmin_ligand_GTP
  )

colnames(chain_hetatm_dis_dimer4)[2:3]<-paste0(colnames(chain_hetatm_dis_dimer4)[2:3],"_dimer4")


anno_dimer4<-merge(interchain_distance_dimer4,chain_hetatm_dis_dimer4,by='Pos')
anno_dimer4[, dimer4_bind := 'no']
anno_dimer4[, dimer4_GTP_bind := 'no']
anno_dimer4[scHAmin_ligand_dimer4 <5 , dimer4_bind := 'yes']
anno_dimer4[GXPMG_scHAmin_ligand_dimer4 <5 , dimer4_GTP_bind := 'yes']



#############################
########dimer4 ##这个是α4-α5/β dimer以B到A的距离
## caculate interchain distance
interchain_distance_dimer4 <- bio3D_minimum_interchain_distances(
  input_file = "7rse.pdb",
  chain_query = "B",
  chain_target = "A")
colnames(interchain_distance_dimer4)[2:3]<-paste0(colnames(interchain_distance_dimer4)[2:3],"_dimer4")

## caculate distance from chain to a specific HETATM
chain_hetatm_dis_dimer4 <- bio3D_chain_hetatm_distances(
  input_file = "6vjj.pdb",
  chain_query = "A",
  hetatm_target = "GNP"
)
colnames(chain_hetatm_dis_dimer4)[2:3]<-paste0(colnames(chain_hetatm_dis_dimer4)[2:3],"_GTP")
chain_hetatm_dis_dimer4 <- chain_hetatm_dis_dimer4 %>%
  rename(
    GXPMG_HAmin_ligand = HAmin_ligand_GTP,
    GXPMG_scHAmin_ligand = scHAmin_ligand_GTP
  )

colnames(chain_hetatm_dis_dimer4)[2:3]<-paste0(colnames(chain_hetatm_dis_dimer4)[2:3],"_dimer4")


anno_dimer4<-merge(interchain_distance_dimer4,chain_hetatm_dis_dimer4,by='Pos')
anno_dimer4[, dimer4_bind := 'no']
anno_dimer4[, dimer4_GTP_bind := 'no']
anno_dimer4[scHAmin_ligand_dimer4 <5 , dimer4_bind := 'yes']
anno_dimer4[GXPMG_scHAmin_ligand_dimer4 <5 , dimer4_GTP_bind := 'yes']





#############################
########dimer3 ##这个是α4-α5/β dimer以A到B的距离
## caculate interchain distance
interchain_distance_dimer3 <- bio3D_minimum_interchain_distances(
  input_file = "7rse.pdb",
  chain_query = "A",
  chain_target = "B")
colnames(interchain_distance_dimer3)[2:3]<-paste0(colnames(interchain_distance_dimer3)[2:3],"_dimer3")

## caculate distance from chain to a specific HETATM
chain_hetatm_dis_dimer3 <- bio3D_chain_hetatm_distances(
  input_file = "6vjj.pdb",
  chain_query = "A",
  hetatm_target = "GNP"
)
colnames(chain_hetatm_dis_dimer3)[2:3]<-paste0(colnames(chain_hetatm_dis_dimer3)[2:3],"_GTP")
chain_hetatm_dis_dimer3 <- chain_hetatm_dis_dimer3 %>%
  rename(
    GXPMG_HAmin_ligand = HAmin_ligand_GTP,
    GXPMG_scHAmin_ligand = scHAmin_ligand_GTP
  )

colnames(chain_hetatm_dis_dimer3)[2:3]<-paste0(colnames(chain_hetatm_dis_dimer3)[2:3],"_dimer3")


anno_dimer3<-merge(interchain_distance_dimer3,chain_hetatm_dis_dimer3,by='Pos')
anno_dimer3[, dimer3_bind := 'no']
anno_dimer3[, dimer3_GTP_bind := 'no']
anno_dimer3[scHAmin_ligand_dimer3 <5 , dimer3_bind := 'yes']
anno_dimer3[GXPMG_scHAmin_ligand_dimer3 <5 , dimer3_GTP_bind := 'yes']



anno_dimer<-merge(anno_dimer1,anno_dimer2,by="Pos",all = T)
anno_dimer<-merge(anno_dimer,anno_dimer3,by="Pos",all = T)
anno_dimer<-merge(anno_dimer,anno_dimer4,by="Pos",all = T)

write.csv(anno_dimer,"anno_dimer.csv")

library(data.table)
anno_final<-fread("anno_final_ma.csv")
anno_final<-merge(anno_dimer,anno_final,by="Pos",all = T)

write.csv(anno_final,"anno_final_all.csv")



















#########################################################
####################PIK3CG and RAF1
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
colnames(chain_hetatm_dis_RAF1)[2:3]<-paste0(colnames(chain_hetatm_dis_RAF1)[2:3],"_GTP")
chain_hetatm_dis_RAF1 <- chain_hetatm_dis_RAF1 %>%
  rename(
    GXPMG_HAmin_ligand = HAmin_ligand_GTP,
    GXPMG_scHAmin_ligand = scHAmin_ligand_GTP
  )

colnames(chain_hetatm_dis_RAF1)[2:3]<-paste0(colnames(chain_hetatm_dis_RAF1)[2:3],"_RAF1")


anno_RAF1<-merge(interchain_distance_RAF1,chain_hetatm_dis_RAF1,by='Pos')
anno_RAF1[, RAF1_bind := 'no']
anno_RAF1[, RAF1_GTP_bind := 'no']
anno_RAF1[scHAmin_ligand_RAF1 <5 , RAF1_bind := 'yes']
anno_RAF1[GXPMG_scHAmin_ligand_RAF1 <5 , RAF1_GTP_bind := 'yes']





## caculate interchain distance
interchain_distance_PIK3CG <- bio3D_minimum_interchain_distances(
  input_file = "1he8.pdb",
  chain_query = "B",
  chain_target = "A")
colnames(interchain_distance_PIK3CG)[2:3]<-paste0(colnames(interchain_distance_PIK3CG)[2:3],"_PIK3")

## caculate distance from chain to a specific HETATM
chain_hetatm_dis_PIK3CG <- bio3D_chain_hetatm_distances(
  input_file = "1he8.pdb",
  chain_query = "B",
  hetatm_target = "GNP"
)
colnames(chain_hetatm_dis_PIK3CG)[2:3]<-paste0(colnames(chain_hetatm_dis_PIK3CG)[2:3],"_GTP")
chain_hetatm_dis_PIK3CG <- chain_hetatm_dis_PIK3CG %>%
  rename(
    GXPMG_HAmin_ligand = HAmin_ligand_GTP,
    GXPMG_scHAmin_ligand = scHAmin_ligand_GTP
  )

colnames(chain_hetatm_dis_PIK3CG)[2:3]<-paste0(colnames(chain_hetatm_dis_PIK3CG)[2:3],"_PIK3")


anno_PIK3CG<-merge(interchain_distance_PIK3CG,chain_hetatm_dis_PIK3CG,by='Pos')
anno_PIK3CG[, PIK3_bind := 'no']
anno_PIK3CG[, PIK3_GTP_bind := 'no']
anno_PIK3CG[scHAmin_ligand_PIK3 <5 , PIK3_bind := 'yes']
anno_PIK3CG[GXPMG_scHAmin_ligand_PIK3 <5 , PIK3_GTP_bind := 'yes']

anno_final<-merge(anno_PIK3CG,anno_RAF1,by="Pos",all = T)


write.csv(anno_final,"anno_final_RAF1_PIK3CG.csv")




