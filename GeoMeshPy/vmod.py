import numpy as np
from math import atan2
from copy import copy, deepcopy
from scipy.spatial import distance
import numpy_indexed as npi
class vertice_modifier:
    
    """
    this class  gets some basic information from the user and returns vertices that are applicable
    for making mesh in GMSH.
    """
    
    def __init__ (self, n_iterations, no_of_faults, all_vers, formations, z_resolution, fault_relation, extent, resolution):
        """
        n_iterations : defines how many geological realizations are going to be used
        no_of_faults : how many faults does exist in the model, if no fault is there set it 0
        all_vers : a list arrays representing all the features, i.e. faults or layers of the model. In case of having
        Fault, mention them first and after them put layers
        formations : an array of formations' names including the last formation (usually called basement in GemPy's terminology)
        z_resolution : for this factor it is needed to have a prediction on the mesh size in the adjacency of layers.
        fault_relation : it defines whether there is a passive fault in the model or not. Refer to Examples to see it
        in more details. See https://github.com/Ali1990dashti/GeoMeshPy
        extent : defines the extent of the model in x, y and z direction.
        resolution : resolution of the model in all direction. If the extent in X direction goes from 0 to 100 m,
        a 20 resolution in x direction says that in every five meter there should be a vertice.
        """
        self.n_iterations = n_iterations # For how many geological realizations you want to make mesh?
        self.no_of_faults = no_of_faults
        self.all_vers = all_vers # vertices of both the layers and faults existing in the model
        self.formations = formations # name of the layers
        self.z_resolution = z_resolution # how fine your mesh will be around contact of layers
        self.fault_relation = fault_relation # relation between faults, are they cutting each other? See
        self.extent = extent
        self.resolution = resolution
        if no_of_faults > 0:
            self.active_F = []
            self.passive_F = []
            for i in range (self.no_of_faults):
                if True not in self.fault_relation[i]:
                    self.passive_F.append(i)
                else:
                    self.active_F.append(i)
        
    @staticmethod
    def chunks(lst, n):
        # A method used for reshaping lists. Almost like numpy.resize().
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    @staticmethod
    def seg_intersect(a1,a2, b1,b2) :
        """
        This method is used when faults are cutting each other. Then, based on this method
        # only important part of the faults will be used for mesh generation.
        """
        da = a2-a1
        db = b2-b1
        dp = a1-b1
        dap = np.empty_like(da)
        dap[0] = -da[1]
        dap[1] = da[0]
        denom = np.dot(dap, db)
        num = np.dot(dap, dp )
        return (num/denom.astype(float))*db + b1

    @staticmethod
    def rotational_sort(list_of_xy_coords, centre_of_rotation_xy_coord, clockwise=True):
        '''
        this function sorts the exterior points in a clockwise order
        '''
        cx, cy = centre_of_rotation_xy_coord
        angles = [atan2(x - cx, y - cy) for x, y in list_of_xy_coords]
        indices = sorted(range(len(angles)), key=angles.__getitem__)
        if clockwise:
            return [list_of_xy_coords[i] for i in indices]
        else:
            return [list_of_xy_coords[i] for i in indices[::-1]]

    def faults_corners(self):
        """
        To include planar fault surfaces in mesh, only four corners of the faults are important because 
        the mesh generator can recreate the planar surface just by its four corners.
        This method runs some calculations on the vertices of the faults and returns four corners and also number of points
        returned for each fault. In cases a fault may need five points rather than four to be recreated correctly.
        """
        faults = self.all_vers[0:self.no_of_faults] # extracts the faults
        faults = [list(column) for column in zip(*faults)] # reshapes the faults 
        # from [[[f1_iter1], [f1_iter2], [f1_iter3]], [[f2_iter1], [f2_iter2], [f2_iter3]]] 
          # to [[[f1_iter1], [f2_iter1]], [[f1_iter2], [f2_iter2]], [[f1_iter2], [f2_iter3]]]
        four_corners = []
        crn = []
        for fault in faults:
            for subfaults in fault:
                nums, counts = np.unique(subfaults[:,1],return_counts=True)
                to_remove_y = subfaults[counts[np.searchsorted(nums, subfaults[:, 1])] < 5]
                remov_ind_y = np.where(np.isin(subfaults,to_remove_y).all(-1)==True)
                b = np.delete(subfaults, (remov_ind_y[0]), axis = 0) # remove redundant points in Y-grid
                nums, counts = np.unique(b[:,-1],return_counts=True)
                to_remove_z = b[counts[np.searchsorted(nums, b[:, -1])] < 5]
                remov_ind_z = np.where(np.isin(b,to_remove_z).all(-1)==True)
                sub_faults = np.delete(b, (remov_ind_z[0]), axis = 0) # remove redundant points in Z-grid
                sorted_sub = sub_faults[np.lexsort((sub_faults[:,1],sub_faults[:,2]))]
                first_slice = sorted_sub[sorted_sub[:,1]==min(sorted_sub[:,1])]
                if np.all(np.diff(first_slice[:,0]) >= 0) or np.all(np.diff(first_slice[:,0]) <= 0):
                    normal = deepcopy(sorted_sub)
                    le_st = sorted_sub[np.where(sorted_sub[:,2] == sorted_sub[0,2])]
                    le_st = len(le_st)
                    le_en = sorted_sub[np.where(sorted_sub[:,2] == sorted_sub[-1,2])]
                    le_en = len(le_en)
                    cor = np.array([sorted_sub[0,:], sorted_sub[int((le_st-1)),:], sorted_sub[-1,:], sorted_sub[-le_en,:]])
                    four_corners.append(cor)
                else:
                    abnormal = deepcopy(sorted_sub)
                    threshold = res_x
                    close_points = abnormal[np.where(np.min(distance.cdist(normal, abnormal),axis=0)<threshold)[0],:]
                    far_po = npi.difference(abnormal, close_points)
                    p_1 = far_po[far_po[:,-1]==max(far_po[:,-1])] # MAX of Z
                    p_2 = p_1[p_1[:,1]==max(far_po[:,1])] # MAX of Y
                    p_3 = p_1[p_1[:,1]==min(far_po[:,1])] # MIN of Y
                    far_repre = np.concatenate([p_3, p_2])

                    sr_far_po = far_po[far_po[:, 1].argsort()]
                    p1 = cor[0][0:3:2]
                    p2 = cor[3][0:3:2]
                    p3 = sr_far_po[0][0:3:2]
                    p4 = sr_far_po[1][0:3:2]
                    fir_p = self.seg_intersect(p1,p2, p3,p4)
                    first_p = np.array([fir_p[0], sr_far_po[0][1], fir_p[1]])
                    p5 = cor[2][0:3:2]
                    p6 = cor[1][0:3:2]
                    p7 = sr_far_po[-1][0:3:2]
                    p8 = sr_far_po[-2][0:3:2]
                    sec_p = self.seg_intersect(p5,p6, p7,p8)
                    second_p = np.array([sec_p[0], sr_far_po[-1][1], sec_p[1]])
                    p_repre = np.array ([first_p, second_p, far_repre[1], far_repre[0]])

                    if np.min (sr_far_po[:,-1]) == np.min (normal[0,-1]): # to handle 5 point faults
                        floor_p = sr_far_po[sr_far_po[:,-1]==min(sr_far_po[:,-1])] # Min of Z: floor points
                        sr_floor = floor_p[floor_p[:, 1].argsort()]
                        if len (sr_floor)>1:
                            p9 = cor[0][0:2]
                            p10 = cor[1][0:2]
                            p11 = sr_floor[0][0:2]
                            p12 = sr_floor[-1][0:2]
                            mid_p = self.seg_intersect(p9,p10, p11,p12)
                            middle_p = np.array([mid_p[0], mid_p[1], sr_floor[0][-1]])
                            p_repre = np.array ([sr_floor[0], middle_p, second_p, far_repre[1], far_repre[0]])
                        else:
                            p_repre = np.array ([cor[0], second_p, far_repre[1], far_repre[0]])
                    four_corners.append(p_repre)
        four_corners=list(self.chunks(four_corners, int (len(four_corners)/self.n_iterations)))
        updated_F = []
        f_c = four_corners
        if len (self.passive_F) > 0:
            for fal_i in range (len(four_corners)):
                for i in range (len(four_corners[fal_i])):
                    if i in self.active_F:
                        updated_F.append(four_corners[fal_i][i])
                    else:
                        a = self.seg_intersect(f_c[fal_i][i][0][0::2], f_c[fal_i][i][3][0::2], f_c[fal_i][0][0][0::2], f_c[fal_i][0][3][0::2])
                        b = self.seg_intersect(f_c[fal_i][i][1][0::2], f_c[fal_i][i][2][0::2], f_c[fal_i][0][1][0::2], f_c[fal_i][0][2][0::2])
                        c = self.seg_intersect(f_c[fal_i][i][1][0::2], f_c[fal_i][i][2][0::2], f_c[fal_i][1][1][0::2], f_c[fal_i][1][2][0::2])
                        d = self.seg_intersect(f_c[fal_i][i][0][0::2], f_c[fal_i][i][3][0::2], f_c[fal_i][1][0][0::2], f_c[fal_i][1][3][0::2])
                        p_a= np.array([[a[0], f_c[fal_i][i][0][1], a[1]],
                                       [b[0], f_c[fal_i][i][1][1], b[1]],
                                       [c[0], f_c[fal_i][i][1][1], c[1]],
                                       [d[0], f_c[fal_i][i][0][1], d[1]]])
                        updated_F.append(p_a)
            updated_F = list(self.chunks(updated_F, int (len(updated_F)/self.n_iterations)))
        else:
            updated_F = four_corners
        lengths_fal = np.array([len(i) for j in updated_F for i in j])
        chang_last = [np.concatenate(d) for d in updated_F]
        len_fal = np.split(lengths_fal, self.n_iterations)
        len_fal = [[i.tolist()] for i in len_fal]
        sub_fourc_list = [l.tolist() for l in chang_last]
        return (sub_fourc_list, len_fal)
        
    
    def contact_generator(self):
        """
        this method is designed to export useable vertices from all the generated ones in modelling
        tools like GemPy. representative_points is another import output of this method.
        This array can help GMSH to assign correct names to each part of your model.
        """
        if self.no_of_faults == 0: # Without any fault, the life is much easier.
            # for such cases you just need to get rid of extra points. See
            # to get an idea about extra points
            cleaned_verti = np.array([])
            vers = [list(column) for column in zip(*self.all_vers)]
            leng = np.array ([])
            representative_points = np.array([])
            for layers in vers:
                for layer_iter in layers:
                    nums, counts = np.unique(layer_iter[:,0],return_counts=True)
                    to_remove_x = layer_iter[counts[np.searchsorted(nums, layer_iter[:, 0])] < 3]
                    remov_ind_x = np.where(np.isin(layer_iter,to_remove_x).all(-1)==True)
                    cleaned_X = np.delete(layer_iter, (remov_ind_x[0]), axis = 0) # remove redundant points in X-grid
                    nums, counts = np.unique(cleaned_X[:,1],return_counts=True)
                    to_remove_y = cleaned_X[counts[np.searchsorted(nums, cleaned_X[:, 1])] < 3]
                    remov_ind_y = np.where(np.isin(cleaned_X,to_remove_y).all(-1)==True)
                    cleaned_Y = np.delete(cleaned_X, (remov_ind_y[0]), axis = 0) # remove redundant points in Y-grid
                    sor = cleaned_Y[np.lexsort((cleaned_Y[:,0],cleaned_Y[:,1]))]
                    leng = np.append (leng, len (sor))
                    cleaned_verti = np.append (cleaned_verti, sor)
                    meds = sor[len(sor)//2]
                    representative_points = np.append(representative_points, meds)
            cleaned_verti = cleaned_verti.reshape(-1,3)
            length_layers = leng.reshape(self.n_iterations,-1)
            cleaned_ver = np.split (cleaned_verti, (np.cumsum(np.cumsum(length_layers,axis=1)[:,-1])).astype('int'))[:-1]
            new_result_list = [l.tolist() for l in cleaned_ver]
            representative_po = representative_points.reshape(self.n_iterations,-1,3)
            representative_p = deepcopy(representative_po)
            representative_p[:,:,-1] += self.z_resolution
            rep_points = np.concatenate((representative_p, representative_p[:,-1:,:]), axis=1)
            rep_points[:,-1,-1] -= (self.z_resolution*2)
            rep_points = np.concatenate(rep_points)
            names = np.tile (self.formations,self.n_iterations).reshape(-1,1)
            repre_pts = np.hstack([rep_points, names]).reshape(self.n_iterations,-1,4).tolist()
            return (new_result_list, length_layers, repre_pts)
        else: # with faults calculations are more intense
            faults =self.all_vers[:self.no_of_faults-len(self.passive_F)]
            faults = [list(column) for column in zip(*faults)]
            sed_layers = self.all_vers[self.no_of_faults:]
            sed_layers = [list(column) for column in zip(*sed_layers)]
            new_result = []
            all_clp = []
            for fault, layers in zip (faults,sed_layers):
                all_fals = np.concatenate(fault[0:])
                resol_x = np.abs(self.extent[1] - self.extent[0]) / self.resolution[0]
                for laye in layers:
                    cp=laye[np.where(np.min(distance.cdist(all_fals, laye),axis=0) < 1.25 * resol_x)[0],:]
                    cleaned_result= npi.difference(laye, cp) # remove undeleted points of layer close to fault
                    # thanks to this step, vertices will be really cut with faults. See 
                    # for more visual clarification in the examples
                    all_clp.append(cp)
                    new_result.append(cleaned_result)
            new_result = list(self.chunks(new_result, int (len(new_result)/self.n_iterations)))
            cleaned_ver = []
            for layers in new_result:
                x_vals = np.array([])
                y_vals = np.array([])
                for sub_layers in layers:
                    nums_x, counts_x = np.unique(sub_layers[:,0],return_counts=True)
                    nums_y, counts_y =np.unique(sub_layers[:,1],return_counts=True)
                    for val_x, freq_x in zip (nums_x, counts_x):
                        if freq_x > max(counts_x)//5:
                            x_vals = np.append(x_vals, val_x)
                    for val_y, freq_y in zip (nums_y, counts_y):
                        if freq_y > max(counts_y)//5:
                            y_vals = np.append(y_vals, val_y)
                    cleand_x = sub_layers[np.in1d(sub_layers[:,0], x_vals)]
                    c = cleand_x[np.in1d(cleand_x[:,1], y_vals)]
                    sor = c[np.lexsort((c[:,0],c[:,1]))]
                    spl_num = np.unique(sor[:,1],return_counts=True)[1]
                    sp_sor = np.split(sor, np.cumsum(spl_num)[:-1])
                    sp_point = []
                    jump = (np.abs(sor[1,-1] - sor[0,-1])+1) * 4
                    for i in range (len(sp_sor)):
                        sp = np.split(sp_sor[i], np.where(np.abs(np.diff(sp_sor[i][:,-1])) > jump)[0]+1)
                        sp_point.append(sp)
                    segs = len (sp_point[0])
                    data_p = np.array([])
                    for i in range (segs):
                        for j in range (len(sp_point)):
                            s = sp_point[j][i]
                            m = s[np.lexsort((s[:,0],s[:,1]))]
                            data_p = np.append(data_p, m)
                    data_p = data_p.reshape(-1,3)
                    data_po = np.split(data_p, np.where(np.diff(data_p[:,1])< 0)[0]+1)
                    cleaned_ver.append([data_po])
            last_sort = cleaned_ver
            reconstructed_arr=[] # this reconstuction step makes some extra points in the borderes of point clouds.
            # We already lost some point in the adjacency of faults. So, they should be reconstructed.
            for realizations in last_sort:
                for layers in realizations:
                    for ind,half_layers in enumerate (layers):
                        last=len(layers)-1
                        coordinates=half_layers[np.lexsort((half_layers[:,0],half_layers[:,1]))]
                        x_grid=coordinates[1,0]-coordinates[0,0]
                        sp_cord = [coordinates[coordinates[:,1]==k] for k in np.unique(coordinates[:,1])]

                        if ind==0: # it gives extra points of FIRST split
                            new_x=[]
                            new_y=[]
                            new_z=np.array([])
                            for i in sp_cord:
                                new_x.append([i[-1,0]+x_grid, i[-1,0]+(2*x_grid), i[-1,0]+(3*x_grid), i[-1,0]+(4*x_grid)])
                                new_y.append(i[0,1])
                                if len (i) >1:
                                    n_z = i[:,-1].copy ()
                                    for _ in range (4):
                                        n = (n_z[-1]-n_z[-2])
                                        n_z = np.append(n_z, n_z[-1]+n)
                                    new_z = np.concatenate ([new_z, n_z[-4:]])   
                                else:
                                    new_z.append(i[-1,-1])
                                    new_z = np.array (np.repeat(new_z, repeats=4, axis=0))
                            new_z = np.array([new_z])
                            new_x = np.array(new_x)
                            new_y = np.array(new_y)
                            new_x_y = []
                            for i in range (len(new_y)):
                                new_x_y.append([new_x[i,0],new_y[i]])
                                new_x_y.append([new_x[i,1],new_y[i]])
                                new_x_y.append([new_x[i,2],new_y[i]])
                                new_x_y.append([new_x[i,3],new_y[i]])
                            new_x_y=np.array(new_x_y)
                            final_first=np.hstack([new_x_y, new_z.T])
                            both_first=np.vstack([coordinates, final_first])
                            all_coordinates_first=both_first[np.lexsort((both_first[:,1],both_first[:,0]))]
                            reconstructed_arr.append(all_coordinates_first)

                        if 0 < ind <last: # it gives extra points of the both sides of the middle splits
                            new_x=[]
                            new_y=[]
                            new_z=np.array([])
                            for ind, i in enumerate (sp_cord):
                                new_x.append([i[-1,0]+x_grid, i[-1,0]+(2*x_grid), i[-1,0]+(3*x_grid), i[-1,0]+(4*x_grid)])
                                new_y.append(i[0,1])
                                if len (i) >1:
                                    n_z = i[:,-1].copy ()
                                    for _ in range (4):
                                        n = (n_z[-1]-n_z[-2])
                                        n_z = np.append(n_z, n_z[-1]+n)
                                    new_z = np.concatenate ([new_z, n_z[-4:]])   
                                else:
                                    ave_depth = 1.5*(np.mean(new_z[0:4])-np.mean(new_z[4:8]))
                                    new_z = np.append (new_z[0:4]+ave_depth, new_z)
            #                         new_z = np.array (np.repeat(np.array(new_z), repeats=4, axis=0))
                            new_z = np.array([new_z])
                            new_x = np.array(new_x)
                            new_y = np.array(new_y)
                            new_x_y = []
                            for m in range (len(new_y)):
                                new_x_y.append([new_x[m,0],new_y[m]])
                                new_x_y.append([new_x[m,1],new_y[m]])
                                new_x_y.append([new_x[m,2],new_y[m]])
                                new_x_y.append([new_x[m,3],new_y[m]])
                            new_x_y = np.array(new_x_y)
                            final_middle_front = np.hstack([new_x_y, new_z.T])
                            both_middle_front = np.vstack([coordinates, final_middle_front]) # Front side

                            new_x = []
                            new_y = []
                            new_z = np.array([])
                            for i in sp_cord:
                                new_x.append([i[0,0]-x_grid, i[0,0]-(2*x_grid), i[0,0]-(3*x_grid), i[0,0]-(4*x_grid)])
                                new_y.append(i[0,1])
                                if len (i) >1:
                                    n_z = i[:,-1].copy ()
                                    for _ in range (4):
                                        n = (n_z[1]-n_z[0])
                                        nm = n_z[0]-n
                                        n_z = np.concatenate(([nm],n_z))
                                    new_z = np.append (n_z[:4], new_z)   
                                else:
                                    ave_depth = 1.5*(np.mean(new_z[0:4])-np.mean(new_z[4:8]))
                                    new_z = np.append (new_z[0:4]+ave_depth, new_z)
            #                         new_z=np.array (np.repeat(np.array(new_z), repeats=4, axis=0))
                            new_z = np.array([new_z])
                            new_x = np.array(new_x)
                            new_y = np.array(new_y)
                            new_x_y = []
                            for i in range (len(new_y)):
                                new_x_y.append([new_x[i,0],new_y[i]])
                                new_x_y.append([new_x[i,1],new_y[i]])
                                new_x_y.append([new_x[i,2],new_y[i]])
                                new_x_y.append([new_x[i,3],new_y[i]])
                            new_x_y = np.array(new_x_y)
                            final_middle_end = np.hstack([new_x_y, new_z.T[::-1]]) # End side  
                            both_middle_end_front = np.vstack([both_middle_front, final_middle_end])
                            all_coordinates_middle = both_middle_end_front[np.lexsort((both_middle_end_front[:,1],both_middle_end_front[:,0]))]    
                            if len (all_coordinates_middle) > 0:
                                reconstructed_arr.append(all_coordinates_middle)                 

                        if ind==last: # it gives extra points of LAST split
                            new_x = []
                            new_y = []
                            new_z = [] 
                            new_z = np.array([])
                            for i in sp_cord:
                                new_x.append([i[0,0]-x_grid, i[0,0]-(2*x_grid), i[0,0]-(3*x_grid), i[0,0]-(4*x_grid)])
                                new_y.append(i[0,1])
                                if len (i) >1:
                                    n_z = i[:,-1].copy ()
                                    for _ in range (4):
                                        n = (n_z[1]-n_z[0])
                                        nm = n_z[0]-n
                                        n_z = np.concatenate(([nm],n_z))
                                    new_z = np.append (n_z[:4], new_z)   
                                else:
                                    new_z.append(i[-1,-1])
                                    new_z=np.array (np.repeat(np.array(new_z), repeats=4, axis=0))
                            new_z = np.array([new_z])                
                            new_x = np.array(new_x)
                            new_y = np.array(new_y)
                            new_x_y = []
                            for m in range (len(new_y)):
                                new_x_y.append([new_x[m,0],new_y[m]])
                                new_x_y.append([new_x[m,1],new_y[m]])
                                new_x_y.append([new_x[m,2],new_y[m]])
                                new_x_y.append([new_x[m,3],new_y[m]])
                            new_x_y = np.array(new_x_y)
                            final_last = np.hstack([new_x_y, new_z.T[::-1]])
                            both_last = np.vstack([coordinates, final_last])
                            all_coordinates_last = both_last[np.lexsort((both_last[:,1],both_last[:,0]))]
                            reconstructed_arr.append(all_coordinates_last)
            sor_arr_rec = []
            meds = [] # this list will keep the middle point of each cloud as a represntative point that can
            # be used by GMSH when you want to map your geological units to volumes created there.
            for i in reconstructed_arr:
                sor_sp=i[np.lexsort((i[:,1].round(decimals=4),i[:,0].round(decimals=4)))].copy()
                sor_arr_rec.append(sor_sp)

                meds.append(sor_sp[int (self.resolution[1, 0]//2 + len(sor_sp)//2)])

            meds = np.split (np.array(meds), self.n_iterations)
            reconstructed_ar = list(self.chunks(sor_arr_rec, int (len(sor_arr_rec)/self.n_iterations)))
            final = [np.concatenate(d) for d in reconstructed_ar]
            new_result_list = [l.tolist() for l in final]
            lengths = np.array([])
            for i in reconstructed_ar:
                for j in i:
                    lengths = np.append(lengths, len(j))
            length_layers = np.split(lengths, self.n_iterations)
            no_lay = len(self.formations)
            mover_fa = self.z_resolution
            repre_pts = []
            for i in meds:
                each_lay = np.split (i, no_lay-1)
                for m in range (len(each_lay)):
                    each_lay[m][:,-1] = each_lay[m][:,-1]+mover_fa
                    if m == len(each_lay)-1:
                        mm = deepcopy(each_lay[m])
                        each_lay.append(mm)
                        each_lay[m+1][:,-1] = each_lay[m+1][:,-1]-(2*mover_fa)
                fin = np.concatenate(each_lay)
                fin = fin.astype('object')
                names = np.repeat (self.formations,len(fin)/no_lay)
                names = names.reshape(len(names),-1)
                arr = np.hstack([fin, names]).tolist()
                repre_pts.append(arr)
            return (new_result_list, length_layers, repre_pts)
