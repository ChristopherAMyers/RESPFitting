import numpy as np
import os
import Elements

class ESPPointGenerator():
    def __init__(self):
        self._BondiRadii = {
            'H': 1.45,
            'C': 1.50,
            'N': 2.00,
            'O': 2.00,
            'P': 2.50,
            'S': 2.20
        }
        self._point_positions = np.array([])

    def _find_neighbors(self, elements, coords, radii):
        ''' Find list of neighbors to each coordinate '''
        nbr_list = []
        for n in range(len(elements)):
            diffs = coords - coords[n]
            norms = np.linalg.norm(diffs, axis=1)
            intersections = np.where(norms <= radii[n] + np.array(radii))[0]
            nbrs = [x for x in intersections if x != n]
            nbr_list.append(nbrs)

        return nbr_list

    def print_xyz(self, out_file_loc):
        with open(out_file_loc, 'w') as file:
            file.write('{:d} \n'.format(len(self._point_positions )))
            file.write('Merz-Kollman ESP points \n')
            for pt in self._point_positions :
                file.write('He {:15.8f}  {:15.8f}  {:15.8f} \n'.format(*pt))

    def gen_MK_points(self, elements, coords, intervals=[], density=5.0, outfile=open(os.devnull, 'w')):
        print("\n Generating ESP Points using Merz-Kollman procedure", file=outfile)

        if len(intervals) == 0:
            intervals = [1.4, 1.6, 1.8, 2.0]
        
        uniq_elms = list(set(elements))
        all_pts = []


        for scale in intervals:
            #   generate sphere points for each element unique type
            scale_point_count = 0
            element_pts = {}
            radii =      [Elements.getRadiiBySymbol(x, 2.0)*scale for x in elements]
            uniq_radii = [Elements.getRadiiBySymbol(x, 2.0)*scale for x in uniq_elms]
            nbr_list = self._find_neighbors(elements, coords, radii)
            for i, radius in enumerate(uniq_radii):
                points = []
                num_theta = np.floor(radius * np.pi * np.sqrt(density))
                dTheta = np.pi/num_theta

                theta_list = np.arange(0.0, np.pi, dTheta)
                sinT = np.sin(theta_list)
                cosT = np.cos(theta_list)

                for t, theta in enumerate(theta_list):
                    if t == 0:
                        points.append(np.array([0, 0, radius]))
                        continue
                    pRad = radius * sinT[t]
                    numPhi = int(np.floor(pRad * 2 * np.pi * np.sqrt(density)))
                    dPhi = 2 * np.pi / numPhi
                    phi = np.arange(0, 2*np.pi, dPhi)[0:numPhi]

                    sphere_pts = np.zeros((numPhi, 3))
                    sphere_pts[:, 0] = sinT[t]*np.cos(phi)*radius
                    sphere_pts[:, 1] = sinT[t]*np.sin(phi)*radius
                    sphere_pts[:, 2] = cosT[t]*radius
                    for pt in sphere_pts:
                        points.append(pt)
                points.append(np.array([0, 0, -radius]))
                element_pts[uniq_elms[i]] = points

            #   apply those spheres to each coordinate
            for n, coord in enumerate(coords):
                points = element_pts[elements[n]] + coord

                #   find points that intersect other spheres
                reject_list = set()
                for nbr in nbr_list[n]:
                    diffs = points - coords[nbr]
                    norms = np.linalg.norm(diffs, axis=1)

                    reject_list = reject_list.union( np.where(norms <= radii[nbr])[0])

                #   remove rejections and add to total list
                keep_points = np.delete(points, list(reject_list), axis=0)
                scale_point_count += len(keep_points)
                for pt in keep_points:
                    all_pts.append(pt)
            print("\t vdW interval {:.2f} includes {:d} ESP points".format(scale, scale_point_count), file=outfile)

        print("\t Total number of ESP points: {:d}\n".format(len(all_pts)), file=outfile)
        self._point_positions = np.array(all_pts)
        return self._point_positions



        
        

    
    

