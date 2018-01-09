#############################################
#  FUNCTION to extract Distance Index V2.1
# this version excludes anonalous intersection objects and prints graphs 
#  over a minimum value of similarity and symmetrical comparison#
#############################
############################################################################################
Dist_Index = function (mapA,mapB,colA,colB,abbreviate=FALSE,plot=FALSE,min_sim=0.5,symmetrical=F, neotropical=F){
  # library(dplyr)
  # library(sf)
  # library(vegan) #just to abbreviate
  # library(maps) #world map if plot=TRUE
  if(plot==TRUE){library(maps);dev.new(); par(mfrow=c(3,7))}
  indexA = which(names(mapA)==colA)
  indexB = which(names(mapB)==colB)
  catsA = sort(unique(mapA[[indexA]]))
  catsB = sort(unique(mapB[[indexB]]))
  nra = length(catsA)
  nrb = length(catsB)
  oldwarnings = getOption("warn")
  options(warn = -1)
  # matrices
  if (abbreviate){
    library(vegan)
    dist_matrix <- matrix( data=NA, ncol=nra, nrow=nrb,
                           dimnames = list(make.cepnames1(catsB,nchar=6), make.cepnames1(catsA,nchar=6)))
  } else  { dist_matrix <- matrix( data=NA, ncol=nra, nrow=nrb, dimnames = list(catsB, catsA))
  }
  
  # loops and individual category maps
  for (w in 1:ifelse(symmetrical,nra,(nra-1))){
    cat_a = mapA[which(mapA[[indexA]] == catsA[w]),]
    # cat_a = mapA[which(mapA[[indexA]] == 'Cercomacroides fuscicauda'),]
    # now doing 'dissolve' using group_by and summarize for every SCINAME fc mapA
    # cat_a = cat_a %>% st_cast(.,'POLYGON')
    cat_a = cat_a %>% group_by(SCINAME)%>%summarize()
    print("################################################")
    print(droplevels(cat_a$SCINAME))
    for (z in ifelse(symmetrical,w+1,1):nrb){
      #mapB -> cat_b
      cat_b = mapB[which(mapB[[indexB]] == catsB[z]),]   # using w+z to fill half matrix, calculating each par just once
 # cat_b = mapB[which(mapB[[indexB]] == 'Cercomacroides nigrescens'),]
      # comparing same spp is useless
      # if (as.character(cat_a[[indexA]]) == as.character(cat_b[[indexB]])) { dist_matrix[z,w] = 1 # not working anymore.. :(
      # next    }
      if (cat_a$SCINAME == cat_b$SCINAME) { dist_matrix[z,w] = 1
      next    }
 
      # # if (class(cat_a[1,]$geometry)[1]=="sfc_MULTIPOLYGON"){}
      
      #'dissolve' using group_by and summarize for mapB
      if (colB=="Name"){cat_b = cat_b %>% group_by(Name)%>%summarize()
      print(droplevels(cat_b$Name))} else  {
        # cat_b = cat_b %>% st_cast(.,'POLYGON')
        cat_b = cat_b %>% group_by(SCINAME)%>%summarize()
        print(paste(droplevels(cat_a$SCINAME),"  vs.  ",droplevels(cat_b$SCINAME)))}
      
      # intersection
      inter = st_intersection(cat_a$geometry,cat_b$geometry)
      # if clause for no intersection
      if (length(inter) == 0){ dist_matrix[z,w] = 0;  next  }
      # if clause for geometry collection class (point+line+polygon !)
      # until figure out how to extract only polygons from within collection object
      if(class(inter)[1] == "sfc_GEOMETRYCOLLECTION") { dist_matrix[z,w] = NA;  next  }
      # {  area_inter = sum(st_area(st_zm(inter[[1]],drop=T))) ; next   }
      # just lines in common, not areas
      if(class(inter)[1] == "sfc_MULTILINESTRING") {dist_matrix[z,w] = NA;  next  } # option is to use zero
      
      # when inter is polygon this works
      area_inter = sum(st_area(st_zm(inter,drop=T)))
      area_a = sum(st_area(st_zm(cat_a,drop=T)))
      area_b = sum(st_area(st_zm(cat_b,drop=T)))
      #writing matrix
      D_calc = (area_inter/area_a)*(area_inter/area_b)
      dist_matrix[z,w] = round(as.numeric(as.character(D_calc)),3)
      print(round(dist_matrix[z,w],3))
      
      #plots
      if ( plot==TRUE){
        if(as.numeric(as.character(D_calc)) > min_sim){
          plot(cat_a,col= rgb(1, 0, 1,0.2),
               xlim=bb_max(cat_a,cat_b)$xlim,ylim=bb_max(cat_a,cat_b)$ylim,
               main=paste(make.cepnames1(cat_a$SCINAME),".", make.cepnames1(cat_b$SCINAME)))
          plot(cat_b,add=T,col=rgb(0, 1, 0,0.2))
          plot(inter,add=T,col=rgb(1, 0, 0,0.2))
          if (neotropical){plot(neotropical$geometry,add=T)} else {map('world',add=T,col=1,lwd=1.5)}
          plot(amaz,lwd=2,add=T)
        }
      }   
    }
  }
  options(warn = oldwarnings)
  (dist_matrix)
}
####################################

###################################################
#####################################   BB_MAX   ##############
bb_max = function (x,y) {
  # to get xlim and ylim just call for obj<-bb_max(A,B)
  #                                   obj$xlim
  #                                   obj$ylim
  
  xy = c('xmin'=min(st_bbox(x)[1],st_bbox(y)[1]),
         'ymin'=min(st_bbox(x)[2],st_bbox(y)[2]),
         'xmax'=max(st_bbox(x)[3],st_bbox(y)[3]),
         'ymax'=max(st_bbox(x)[4],st_bbox(y)[4]))
  
  return (list (xy = xy, xlim = c(xy[1],xy[3]), ylim = c(xy[2],xy[4]) ))
  
}
###################################################################################################

#######################################################################################################################################
#######################################################################################################################################
#       CO-OCCURRENCE PATTERNS FUNCTION
#######################################################################################################################################
#######################################################################################################################################
co_occurrence_patterns = function (D_df=D_df,min_sim=0.7, method = 'vector_sim',
              map=map,plot_map=TRUE, ratio_IU=0.5, min_sim_vect = 0.7,
              min_sim_to_pattern = 0.5, save_spat_obj=FALSE, vector_sim=vector_sim,
              spatial_check = TRUE) {
  
####### check spatial function  ############# #####
check_spatial = function(...) {   # args = map_sp, map_sp_2, map_sp_pattern, min-sim - already running in co-occ function
  # areas of polygons
      asp = st_area(st_geometry(map_sp))
      asp_2 = st_area(st_geometry(map_sp_2))
      asp_pattern = st_area(st_geometry(map_sp_pattern))
      ai = st_area(st_intersection(st_geometry(map_sp),st_geometry(map_sp_2)))
      ai_pattern_sp = st_area(st_intersection(st_geometry(map_sp),st_geometry(map_sp_pattern)))
      ai_pattern_sp_2 = st_area(st_intersection(st_geometry(map_sp_2),st_geometry(map_sp_pattern)))
  # D_similarities    
      D_chk =  (as.numeric((ai/asp)*(ai/asp_2)) < min_sim)  # suposely to be equall to D_loop
      D_chk_pattern = (as.numeric((ai_pattern_sp/asp)*(ai_pattern_sp/asp_pattern))<min_sim |
                         as.numeric((ai_pattern_sp_2/asp_2)*(ai_pattern_sp_2/asp_pattern))<min_sim)
   # test   
      if (D_chk | D_chk_pattern) { 
        if (D_chk & D_chk_pattern) { check_spat = 'D_chk < min_sim & D to pattern sp < min-sim'}
        if (D_chk) { check_spat = 'D_chk < min_sim'}
        if (D_chk_pattern)  { check_spat = c('D to pattern sp < min-sim') }
      }   else     {   check_spat = TRUE   }
      return (check_spat)
}
    
##### library, decription, warnings etc ##### #######
  # D_df is the simetrical data.frame (matrix) of similarities between species
  # sim_min is the threshold similarity value to comparison between species
  # method 'standard' uses a fixed threshold; 'adaptative' makes an equalization based on species' areas of occurrence
   # list of results creation!!
    
    library(dplyr)
    library(sf)
    oldw <- getOption("warn")
    options(warn = -1)

#### output settings & meta-info ############ ####
    output = list() # output
    output[['meta_info']] = paste('Species co-occurrence patterns based on geographical distribution similarity', '\n',
                                'Presets:','\n',
                                'minimal similarity = ', min_sim, '\n',
                                'min similarity to build vector of similarity between species = ',min_sim_vect, '\n',
                                'min similarity to first(pattern) species = ', min_sim_to_pattern, '\n',
                                'min ratio between cumulative intersection and union = ', ratio_IU, '\n'
                                )
############ FILTER etc. #################### #############
    if (any(colnames(D_df) != rownames(D_df))){stop("Error: matrix
                            not symmetrical|rows and columns are not identical")}
    for(z in 1:nrow(D_df)){D_df[z,z] = NA} # transform all 1 in NA for equal species in row/col
    D_df = ifelse(D_df > min_sim,D_df,NA)   # all D_df < min_sim are converted to NA
    filt_r_min = apply(D_df,1,function(x) all(x<min_sim))  # all columns (spp similarities) < min_sim?
    D_df = D_df[-is.na(filt_r_min),-is.na(filt_r_min)]   # only spp with at least one similarity > min
    nome = rownames(D_df)
    ratio_inter_union = 1
##### vector_D_sim ########################## ####
    #only with means of non-NA similarity above min_sim_vect threshold
    if (method == 'vector_sim'){vector_D_sim = vector_sim} else {
      vector_D_sim = order(rowMeans(ifelse(D_df > min_sim, D_df, NA), na.rm = TRUE), decreasing = TRUE)  # ordered vector with D_sim means by sp
    }
    if (polyg) {  x11();  par(mfrow = c(2,3), mar = c(1,3,1,2))  }
        
############ FOR LOOP ###########################
for (sp in vector_D_sim) {
 # n=1
 # sp = vector_D_sim[n]; n=n+1
####### first settings ###################### ##########
        # sp_col = which(!is.na(D_df[sp,]), arr.ind = T)  # can be which > min_sim
        sp_col = which(D_df[sp,] > min_sim, arr.ind = T)
        vector_Dsp = D_df[sp,][sp_col]  # these next two lines are suposed to be outside the internal loop (while D > min_sim)...think about
        sp_2 = which (D_df[sp,] == max(D_df[sp,],na.rm=T))[1] # required if there are more than one spp with max D_sim
        D_loop = D_df[sp,sp_2]
        sp_pattern = sp # keeps the #1 pattern species index
        pattern_out = paste('pattern_',nome[sp_pattern],sep="")
        output[[pattern_out]] = list()
        loop = 'go'
        n_loop = 0
        nxt = 0
        ratio_inter_union = 1
        cat(pattern_out," Starts",'\n')
      
        inter_comb = st_sf(id = 1, geometry = st_sfc(lapply(1, function(x) st_multipolygon())))
        union_comb = st_sf(id = 1, geometry = st_sfc(lapply(1, function(x) st_multipolygon())))
##########  WHILE / REPEAT LOOP #################
  repeat { 
####### first BREAK STATEMENT  ############## ########
# estou com d?vida cruel:
# se sp for escolhida de novo, ap?s algumas rodadas, como sp, atrav?s de outra esp?cie
    # que nao tenha sido comparada ainda com sp... ser? que isto est? previsto nas condi??es
    # de break?
    if ((D_df[sp,sp_2] < min_sim) | is.na(D_df[sp,sp_2]) ) { break; print('1st break assessment - D[sp,sp_2]<min_sim') } # meaning there is no sp_2 with signficant D_sim to sp_pattern; folows the next species in for loop 
# migrate this part to the end
    # if ((loop == 'break') | (D_loop < min_sim) | (ratio_inter_union < ratio_IU)) { break } 
    # if ( sp != sp_pattern) { if ( (D_df[sp_pattern,sp] < min_sim_to_pattern) |
    #    (D_df[sp_pattern,sp_2] < min_sim_to_pattern) |   # Dist to original species must be above min_sim_to_pattern
    #    (is.na(D_df[sp_pattern,sp])) |
    #    (is.na(D_df[sp_pattern,sp_2]))){ break }}
    
########## Setup vector D_sp / D_sp_2  ###### ######
        n_loop = n_loop + 1
        sp_col_2 = which(!is.na(D_df[sp_2,]), arr.ind = T); #sp_col_2
        vector_Dsp_2 = D_df[sp_2,][sp_col_2]; #vector_Dsp_2
# if sp changes after first loop, recalculates vector_Dsp
        if (sp != sp_pattern) { sp_col = which(!is.na(D_df[sp,]), arr.ind = T)  
                                                vector_Dsp = D_df[sp,][sp_col]      }
#########  Data Frame  Update  ############## #########
        # out_table = 'unchanged'
    #  comparar apenas com primeira coluna 'spec', tanto sp quanto sp_2
    if (!(sp %in% output[[pattern_out]][['df']][,'spec_ind']) |
        !(sp_2 %in% output[[pattern_out]][['df']][,'spec_ind'])) { # sp & sp_2 %in% spec
      
      if ((sp %in% output[[pattern_out]][['df']][,'spec_ind']) &
          !(sp_2 %in% output[[pattern_out]][['df']][,'spec_ind'])) {  # sp in spec but sp_2 not
          spec = c(rep(nome[sp_2],length(vector_Dsp_2)))
          spec_ind = c(rep(sp_2,length(vector_Dsp_2)))
          spec_2 = c(nome[sp_col_2])
          spec_ind_2 = sp_col_2
          D_sim = as.numeric(vector_Dsp_2)
          used = numeric(length = length(sp_col_2)) 
          # D_sp_to_patt = ifelse(sp != sp_pattern, D_df[sp,sp_pattern],1)
          # D_sp_2_to_patt = ifelse(sp_2 != sp_pattern, D_df[sp_2,sp_pattern],1)
          }
      
      if (!(sp %in% output[[pattern_out]][['df']][,'spec_ind']) &
          (sp_2 %in% output[[pattern_out]][['df']][,'spec_ind'])) { # sp_2 in spec but not sp
          spec = c(rep(nome[sp],length(vector_Dsp)))
          spec_ind = c(rep(sp,length(vector_Dsp)))
          spec_2 = nome[sp_col]
          spec_ind_2 = sp_col
          D_sim = as.numeric(vector_Dsp)
          used = numeric(length = length(sp_col))
          # D_sp_2_to_patt = ifelse(sp_2 != sp_pattern, D_df[sp_2,sp_pattern],1)  
          }
      
      if (!(sp %in% output[[pattern_out]][['df']][,'spec_ind']) &
          !(sp_2 %in% output[[pattern_out]][['df']][,'spec_ind'])) { 
        
          spec = c(rep(nome[sp],length(vector_Dsp)),rep(nome[sp_2],length(vector_Dsp_2)))
          spec_ind = c(rep(sp,length(vector_Dsp)),rep(sp_2,length(vector_Dsp_2)))
          spec_2 = c(nome[sp_col],nome[sp_col_2]) # B.T.W == c(names(vector_Dsp),names(vector_Dsp_2)))
          spec_ind_2 = c(sp_col,sp_col_2)
          D_sim = as.numeric(c(vector_Dsp, vector_Dsp_2))
          used = numeric(length = length(sp_col)+length(sp_col_2))
          # D_sp_to_patt = ifelse(sp != sp_pattern, D_df[sp,sp_pattern],1)      
          }
      
      tab = data.frame(spec,spec_2,spec_ind,spec_ind_2,D_sim = as.numeric(D_sim),used=as.numeric(used))#, D_sp_to_patt=as.numeric(D_sp_to_patt),D_sp_2_to_patt=as.numeric(D_sp_2_to_patt))
      rownames(tab)= NULL
      output[[pattern_out]][['df']] = unique(rbind(as.data.frame(output[[pattern_out]][['df']]),tab)) 
      # out_table = 'changed'  # forgot why it is needed
    }
         
##########   CALCULUS  POLYGONS ###########   #################
      if (n_loop ==1) {ap = which(tab_abrev$sp_abrev == nome[sp_pattern])
      map_sp_pattern = map[ap,"SCINAME"] %>% st_cast("MULTIPOLYGON")
      map_sp_pattern = map_sp_pattern %>% group_by('SCINAME') %>% summarise()
      }
      ab = which(tab_abrev$sp_abrev == nome[sp])
      map_sp = map[ab,"SCINAME"] %>% st_cast("MULTIPOLYGON")
      map_sp = map_sp %>% group_by('SCINAME') %>% summarise()# %>% st_cast("POLYGON")
      ab_2 = which(tab_abrev$sp_abrev == nome[sp_2])
      map_sp_2 = map[ab_2,"SCINAME"] %>% st_cast("MULTIPOLYGON")
      map_sp_2 = map_sp_2 %>% group_by('SCINAME') %>% summarise()
      if (!st_is_valid(map_sp_pattern)) {map_sp_pattern = st_buffer(map_sp,0)}
      if (!st_is_valid(map_sp)) {map_sp = st_buffer(map_sp,0)}
      if (!st_is_valid(map_sp_2)) {map_sp_2 = st_buffer(map_sp_2,0)}
      # map_sp = st_make_valid(map_sp) # not working using workaround 'buffer(x,0)' (needs lwgeos and this package is not working with sf)
      # map_sp_2 = st_make_valid(map_sp_2)
      inter_sp = st_intersection(st_geometry(map_sp), st_geometry(map_sp_2))
      if (!st_is_valid(inter_sp)) {inter_sp = st_buffer(inter_sp,0)}
      union_sp = st_union(st_geometry(map_sp), st_geometry(map_sp_2))
      if (!st_is_valid(union_sp)) {union_sp = st_buffer(union_sp,0)}
      st_crs(union_comb) = st_crs(inter_comb) = st_crs(map_sp)
      # combine inter and union
      if (n_loop == 1) {inter_comb = inter_sp; union_comb = union_sp} else {
        inter_comb = st_intersection(st_geometry(inter_comb), st_geometry(inter_sp))
        union_comb = st_union(st_geometry(union_comb), st_geometry(union_sp))
        # checkin and fixing polygon validity
        if (!st_is_valid(inter_comb)) {inter_comb = st_buffer(inter_comb,0)}
        if (!st_is_valid(union_comb)) {union_comb = st_buffer(union_comb,0)} }
      # ratio used to break pattern
      ratio_inter_union = as.numeric(as.character(st_area(inter_comb)/st_area(union_comb)))
      # leave the loop condition
      if (ratio_inter_union < min_sim) {break; print('Calculus polygons: ratio_inter_union < min_sim')}

######### plotS ############################# ###################
      if (plot_map){
        if(n_loop ==1){ 
        x = c(st_bbox(union_comb)[1] - 7 , st_bbox(union_comb)[3] + 7)  # aumentar a 'folga' aqui no plot
        y = c(st_bbox(union_comb)[2] - 5 , st_bbox(union_comb)[4] + 5)
        plot(union_sp,col=rgb(0,.9,0,.2), xlim = x, ylim = y,
             main=paste(pattern_out))#,'ratio =',round(ratio_inter_union,2)))
        plot(map_sp,col=rgb(.3,.4,0,.2),add=T,lwd=3,lty=5)
        plot(map_sp_2,col=rgb(.1,.2,.8,.2),add=T)
        plot(inter_sp,col=rgb(1,0,0,.2),add=T)
        plot(neotropical$geometry,add=T,lwd=1.5)
        factor_correc_pos = 1.3 - (n_loop * 0.1)
        text(x=(x[1]+ 0.6*(x[2]-x[1])), y= (y[1] + 1.3*(y[2]-y[1])),
             paste(nome[sp_2],'D =',D_df[sp,sp_2]),col='blue')#,adj=c(1,1))
        } else  {
        plot(union_comb, col= rgb(0,1,0,.3),add=T,lwd=2, main="testeset")
        plot(map_sp_2,col=rgb(.1,.2,.8,.2),add=T,lty=n_loop)
        plot(inter_comb, col=rgb(1,0,0,.2),add=T,lwd=2)
        factor_correc_pos = 1.3 - (n_loop * 0.1)
        if (sp != sp_pattern){text(x=(x[1]+ 0.6*(x[2]-x[1])), y= (y[1] + factor_correc_pos*(y[2]-y[1])), 
                                   paste(nome[sp],"/",nome[sp_2],'D =',D_df[sp,sp_2]),col='blue')}#,adj=c(1,1))}
        # text(x=(x[1]+ 0.6*(x[2]-x[1])), y= (y[1] + factor_correc_pos*(y[2]-y[1])), paste(nome[sp],nome[sp_2],'D =',D_loop),col='red',adj=c(0,.5))
        }
      }
######## TAB SEQUENCE comparisons  - output   ############################ 
      print(paste(nome[sp],' vs. ', nome[sp_2], ': ',D_df[sp,sp_2],sep=''))
      
      tab_seq = data.frame(rbind(data.frame(spec = nome[sp], spec_2 = nome [sp_2], spec_ind = sp,
          spec_ind_2 = sp_2, D_sim = D_df[sp,sp_2], ratio_IU_cum = ratio_inter_union),data.frame(spec = nome[sp_2], spec_2 = nome [sp],
          spec_ind = sp_2, spec_ind_2 = sp, D_sim = D_df[sp,sp_2],ratio_IU_cum = ratio_inter_union)))
      rownames(tab_seq)= NULL
      
      if (n_loop == 1) {output[[pattern_out]][['comparison_sequence']] = tab_seq} else {
      output[[pattern_out]][['comparison_sequence']] =
        data.frame(unique(rbind(output[[pattern_out]][['comparison_sequence']],tab_seq)))
      }
      
      # mark already analysed
      current_pair_index =  with( output[[pattern_out]][['df']], 
                                  which( (spec_ind == sp & spec_ind_2 == sp_2) | (spec_ind == sp_2 & spec_ind_2 == sp)))# & output[[pattern_out]][['df']]$spec_ind_2 == sp_2)]
      
      output[[pattern_out]][['df']][as.numeric(current_pair_index),'used'] = n_loop
      
######### NEXT SPP  CHOICE ################# ################
      # copy dataframes
      d = output[[pattern_out]][['df']][with(output[[pattern_out]][['df']],which(D_sim > min_sim_vect & used == 0)),] # only unused pairs are selected
      s = output[[pattern_out]][['comparison_sequence']] # already compared for check
      D_table_loop = unique(d[order(d[,'D_sim'],decreasing = T),])
      nxt = 1
      if (nrow (D_table_loop) < 1) {break; print('next spp choice: no more combinations above min_sim')}
      new_sp = D_table_loop$spec_ind[nxt]
      new_sp_2 = D_table_loop$spec_ind_2[nxt]
      if (is.na(new_sp) | is.na(new_sp_2)) {break; print('next spp choice: is NA new_sp or new_sp_2')}
      cat(nome[new_sp],new_sp,' X ', nome[new_sp_2], new_sp_2, '; D =', D_df[new_sp,new_sp_2],'\n')
########## loop for repeated new pairs  ##### #######
      while (any( new_sp == s$spec_ind & new_sp_2 == s$spec_ind_2 ) | 
                 any( new_sp_2 == s$spec_ind & new_sp == s$spec_ind_2 ) ){   #in case this pair was already analysed
        
        repeated_pair_index =  with( output[[pattern_out]][['df']], 
                                    which( (spec_ind == new_sp & spec_ind_2 == new_sp_2) | (spec_ind == new_sp_2 & spec_ind_2 == new_sp)))# & output[[pattern_out]][['df']]$spec_ind_2 == sp_2)]
        output[[pattern_out]][['df']][as.numeric(repeated_pair_index),'used'] = n_loop
        nxt = nxt + 1
        if (nxt > length (D_table_loop)) {loop = 'break';break; print('loop for repeated new pairs: no more combinations above min_sim')}
        new_sp = D_table_loop$spec_ind[nxt]
        new_sp_2 = D_table_loop$spec_ind_2[nxt]
        paste(nome[new_sp],new_sp,' \ ', nome[new_sp_2], new_sp_2, 'D=', D_df[new_sp,new_sp_2])
      }
######## new_sp / new_sp_2 to sp / sp_2 ##### ####
      if (loop == 'break') {  break ; print('loop for repeated new pairs/new_sp to sp: no more combinations above min_sim') }  else  {
                
        sp = new_sp; sp_2 = new_sp_2; D_loop = D_df[sp,sp_2]    }
      
      # output;sp;sp_2
      paste(nome[new_sp],new_sp,' X ', nome[new_sp_2], new_sp_2, '; D =', D_df[new_sp,new_sp_2], '/ loop =', n_loop)
      
####### second BREAK AsSEssMENT  ############ ########
      # just to check loop manualy # check loop conditions
      if ((D_df[sp,sp_2] < min_sim) | is.na(D_df[sp,sp_2]) ) { break; cat('second BREAK AsSEsMENT: D_loop sp/sp_2 below min_sim') } # meaning there is no sp_2 with signficant D_sim to sp_pattern; folows the next species in for loop
      if ((loop == 'break') | (D_loop < min_sim) | (ratio_inter_union < ratio_IU)) {break; cat('second BREAK AsSEsMENT: (loop == break) | (D_loop < min_sim) | (ratio_inter_union < ratio_IU) ') }
      if ( sp != sp_pattern) { 
        if ( (D_df[sp_pattern,sp] < min_sim_to_pattern) |
         (D_df[sp_pattern,sp_2] < min_sim_to_pattern) |   # Dist to original species must be above min_sim_to_pattern
         (is.na(D_df[sp_pattern,sp])) |
         (is.na(D_df[sp_pattern,sp_2]))) { 
          if (spatial_check) {
          if (check_spatial()) {print('second BREAK AsSEsMENT: spatial test OK ') } else {  break; print('D_sim with bad value but spatial test D<min_sim ' ) }     #if there is problem with dist to pattern, check must be done spatially...
          }   }
         }
# end loop repeat   ####
    }
######   SAVING union / intersection ######## ########
    if (save_spat_obj){
      if (!is.na(inter_comb))  output[[pattern_out]][['intersection']]=inter_comb
      if (!is.na(union_comb))  output[[pattern_out]][['union']] = union_comb   
        }
  
    cat(pattern_out," Ends",'\n') 
    output[[pattern_out]][['comparison_sequence']]
    
# end loop for  #####
  }
  options(warn = oldw)
return(output)
# end function
} 

####################################################################
# Filter maps by overlay with vector
###################################################################
filter_by_overlay = function(mapA,mask){  # removed categories vector
  catsA = unique(as.character(mapA$SCINAME))
  # overlay_matrix <- matrix( data=NA, ncol=3, nrow=length(catsA),
  #                           dimnames = list(catsA,c("species", "amazon","prop_amazon")))
  overlay_df = data.frame('species'=NA,'prop_amaz'=NA,'area_sp'=NA)
  for (n in 1:length(catsA)) {  # 216 unique spp. from 'map'$SCINAME
    index = which(mapA$SCINAME == catsA[n])
    map_sp = mapA[index,] # column where the categories are specified
    # map_sp = furna_map_amaz[which(furna_map_amaz$SCINAME == furna_map_amaz_spp[n]),]
    map_sp = map_sp %>% group_by(SCINAME)%>%summarize()
    # print(as.character(droplevels(map_sp$SCINAME[n])))
    # overlay_matrix[n,1] = as.character(catsA[n])#as.character(droplevels(map_sp$SCINAME[n]))
    overlay_df[n,1]=as.character(droplevels(map_sp$SCINAME))
    #intersection
    inter = st_intersection(st_geometry(map_sp),st_geometry(mask))
    # if clause for no intersection
    # if (length(inter) == 0){ overlay_matrix[n,2:3] = NA;  next  }
    # if(class(inter)[1] == "sfc_GEOMETRYCOLLECTION") {overlay_matrix[n,2:3] = NA;  next  }
    # if(class(inter)[1] == "sfc_MULTILINESTRING") {overlay_matrix[n,2:3] = NA;  next } 
    
    if (length(inter) == 0){ overlay_df[n,2] = NA;  next  }
    if(class(inter)[1] == "sfc_GEOMETRYCOLLECTION") {overlay_df[n,2] = NA;  next  }
    if(class(inter)[1] == "sfc_MULTILINESTRING") {overlay_df[n,2] = NA;  next } 
    
    # # when inter is polygon this works
    area_inter = sum(st_area(st_zm(inter,drop=T))); area_inter = sum(area_inter)
    area_sp = sum(st_area(st_zm(map_sp,drop=T))); area_sp = sum(area_sp)
    overlay_df[n,3] = round(area_sp ,4)
    #writing matrix
    overlay = as.numeric(area_inter/area_sp)
    # overlay_matrix[n,3] = as.numeric(round(overlay,4))
    if (is.numeric(overlay)){overlay_df[n,2] = round(overlay,4)} else {overlay_df[n,2] = NA}
    print(paste("[",n,"] ",overlay_df[n,1],": Overlay_amazon >",prop,"/ = ", overlay_df[n,2]))
  }
  return(data.frame(overlay_df))
}
####################################################################

########################
### function adapted from vegan 'make.cepnames' - changed to use _ between abbreviated names...
make.cepnames1=function (names, nchar=4, seconditem = FALSE) {
  names <- make.names(names, unique = FALSE)
  names <- gsub("\\.[\\.]+", ".", names)
  names <- gsub("\\.$", "", names)
  names <- lapply(strsplit(names, "\\."), function(x) if (length(x) > 
                                                          1) 
    substring(x, 1, nchar)
    else x)
  names <- unlist(lapply(names, function(x) if (length(x) > 
                                                1) 
    paste(x[c(1, if (seconditem) 2 else length(x))], collapse = "_")
    else x))
  names <- abbreviate(names, 2*nchar+1)
  names <- make.names(names, unique = TRUE)
  names
}

###############################
#############################################################
# create a "map" sf object with all features of interest
########### MAKE MAP FUNCTION    ############################################
make_map = function (list_names,dsn){
  for (i in 1:length(list_names)){
    if (i == 1) {map = st_read(dsn = dsn,layer = list_names[i])}
    else {map = rbind(map, st_read(dsn=dsn,layer=list_names[i]))}
  }
  return(map)
}
#######################

################## ggplot function wrap   ##########
ggp = function(data,column){
  x11()
  ggplot()+geom_sf(data=data,aes(fill=factor(column)))
}
###########################################################



#### MAKING MAPS ####

map_family = function (family, dir = "C:/SIG2015/BIRDLIFE_10/",functions_in=TRUE )
# libraries
library(sf)
##### auxiliary functions ##########################################
if (functions_in) {
make_map = function (list_names,dsn){
  for (i in 1:length(list_names)){
    if (i == 1) {map = st_read(dsn = dsn,layer = list_names[i])}
    else {map = rbind(map, st_read(dsn=dsn,layer=list_names[i]))}
  }
  return(map)
}
####################################################################
filter_by_overlay = function(mapA,mask){  # removed categories vector
  catsA = unique(as.character(mapA$SCINAME))
  # overlay_matrix <- matrix( data=NA, ncol=3, nrow=length(catsA),
  #                           dimnames = list(catsA,c("species", "amazon","prop_amazon")))
  overlay_df = data.frame('species'=NA,'prop_amaz'=NA,'area_sp'=NA)
  for (n in 1:length(catsA)) {  # 216 unique spp. from 'map'$SCINAME
    index = which(mapA$SCINAME == catsA[n])
    map_sp = mapA[index,] # column where the categories are specified
    # map_sp = furna_map_amaz[which(furna_map_amaz$SCINAME == furna_map_amaz_spp[n]),]
    map_sp = map_sp %>% group_by(SCINAME)%>%summarize()
    # print(as.character(droplevels(map_sp$SCINAME[n])))
    # overlay_matrix[n,1] = as.character(catsA[n])#as.character(droplevels(map_sp$SCINAME[n]))
    overlay_df[n,1]=as.character(droplevels(map_sp$SCINAME))
    #intersection
    inter = st_intersection(st_geometry(map_sp),st_geometry(mask))
    # if clause for no intersection
    # if (length(inter) == 0){ overlay_matrix[n,2:3] = NA;  next  }
    # if(class(inter)[1] == "sfc_GEOMETRYCOLLECTION") {overlay_matrix[n,2:3] = NA;  next  }
    # if(class(inter)[1] == "sfc_MULTILINESTRING") {overlay_matrix[n,2:3] = NA;  next } 
    
    if (length(inter) == 0){ overlay_df[n,2] = NA;  next  }
    if(class(inter)[1] == "sfc_GEOMETRYCOLLECTION") {overlay_df[n,2] = NA;  next  }
    if(class(inter)[1] == "sfc_MULTILINESTRING") {overlay_df[n,2] = NA;  next } 
    
    # # when inter is polygon this works
    area_inter = sum(st_area(st_zm(inter,drop=T))); area_inter = sum(area_inter)
    area_sp = sum(st_area(st_zm(map_sp,drop=T))); area_sp = sum(area_sp)
    overlay_df[n,3] = round(area_sp ,4)
    #writing matrix
    overlay = as.numeric(area_inter/area_sp)
    # overlay_matrix[n,3] = as.numeric(round(overlay,4))
    if (is.numeric(overlay)){overlay_df[n,2] = round(overlay,4)} else {overlay_df[n,2] = NA}
    print(paste("[",n,"] ",overlay_df[n,1],": Overlay_amazon >",prop,"/ = ", overlay_df[n,2]))
  }
  return(data.frame(overlay_df))
}
####################################################################
}


# BirdLife 9.1 species database
BL = read.csv('C:/SIG2018/BIRDLIFE_10/BirdLife_v9.1.csv', header=T)

# generical maps
neotropical = st_read(dsn='C:/SIG2018/CO_OCCUR_PATT/shapefiles_amazon',layer="Lowenberg_Neto_2014")
amazonia = st_read(dsn='C:/SIG2018/CO_OCCUR_PATT/shapefiles_amazon',layer= "amazonia_polygons")  ### need to copy this shape !!!!!
amaz = st_read(dsn='C:/SIG2018/CO_OCCUR_PATT/shapefiles_amazon',layer="amazonia_sensulatissimo")
shapefiles_folder = 'C:/SIG2018/BIRDLIFE_10/'
family_dir = paste(shapefiles_folder,family,sep='')

########################################################################
# modified mask to amazon excluding 'Andes' (> ~ 650m)
a = which(amazonia$name == "Andes" )
am = amazonia[-a,]
am = am%>% group_by(name) %>% summarise()  # plot(st_geometry(am),lwd=3,col=rgb(1,0,0,.3), add=T)

# FAMILY list shp files BL10
fur_shpBL10 = list.files(family_dir, pattern = ".shp")
fur_shpBL10 = fur_shpBL10[- grep(".shp.xml", fur_shpBL10)]
fur_shpBL10 = gsub( ".shp", "", fur_shpBL10)
####   FAMILIES MAPS   ####
# Lists and selection of species with intersection with Amazon 
furna_map = make_map (list_names = fur_shpBL10,dsn=family_dir)  # function make_map above
st_crs(amaz)=st_crs(furna_map)
furna_map_amaz = furna_map[amaz,] # just SCINAME column
furna_map_spp = unique(furna_map$SCINAME)
furna_map_amaz_spp = unique(furna_map_amaz$SCINAME)
furna_map_amaz_spp=droplevels(furna_map_amaz_spp)  # how many spp 'touching' amazon sensu lato

# Family members in BL9.1 table
fur_BL9 = BL[which(BL$Family == family),]  # LIST OF taxa in BL9.1 check-list
# Family taxons not in BL9.1 or not in shapefiles (don't know exactly)
tax_BL9_not_10 = setdiff(as.vector(fur_BL9$SCINAME),furna_map_spp) # usefull for future check

#filter
furna_filter_amaz = filter_by_overlay(mapA=furna_map_amaz,mask=am)
names(furna_filter_amaz)
# filter 0.75 
filter.75 = which(furna_filter_amaz$prop_amaz > .75)
furna_map_.75_spp = furna_filter_amaz[filter.75,]
nrow(furna_map_.75_spp)  # 110 Furnariidae spp > 0.75 % distrib in amazon sensu lato
# filter 0.6
filter.6 = which(furna_filter_amaz$prop_amaz > .6)
furna_map_.6_spp = furna_filter_amaz[filter.6,]
nrow(furna_map_.6_spp)  # 117 Furnariidae spp > 0.6 % distrib in amazon sensu lato
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
furna_map_.75 = furna_map_amaz[furna_map_amaz$SCINAME %in% furna_map_.75_spp[,1],]
furna_map_.6 = furna_map_amaz[furna_map_amaz$SCINAME %in% furna_map_.6_spp[,1],]
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
##


