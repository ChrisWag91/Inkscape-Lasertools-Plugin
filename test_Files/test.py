def point_inside_csp(p, csp, on_the_path=True):

    x, y = p
    ray_intersections_count = 0

    for subpath in csp:
        for i in range(1, len(subpath)):
            sp1, sp2 = subpath[i-1], subpath[i]
            ax, bx, cx, dx = csp_parameterize(sp1, sp2)[::2]
            if ax == 0 and bx == 0 and cx == 0 and dx == x:
                # we've got a special case here
                b = csp_true_bounds([[sp1, sp2]])
                if b[1][1] <= y <= b[3][1]:
                    # points is on the path
                    return on_the_path
                else:
                    # we can skip this segment because it wont influence the answer.
                    pass
            else:
                for t in csp_line_intersection([x, y], [x, y+5], sp1, sp2):

                    if t == 0 or t == 1:
                        # we've got another special case here
                        y1 = csp_at_t(sp1, sp2, t)[1]

                        if y1 == y:
                            # the point is on the path
                            return on_the_path
                        # if t == 0 we sould have considered this case previously.

                        if t == 1:
                            # we have to check the next segmant if it is on the same side of the ray
                            st_d = csp_normalized_slope(sp1, sp2, 1)[0]

                            if st_d == 0:
                                st_d = csp_normalized_slope(sp1, sp2, 0.99)[0]

                            for j in range(1, len(subpath)+1):

                                if (i+j) % len(subpath) == 0:
                                    continue  # skip the closing segment

                                sp11, sp22 = subpath[(
                                    i-1+j) % len(subpath)], subpath[(i+j) % len(subpath)]
                                ax1, bx1, cx1, dx1 = csp_parameterize(
                                    sp1, sp2)[::2]

                                if ax1 == 0 and bx1 == 0 and cx1 == 0 and dx1 == x:
                                    continue  # this segment parallel to the ray, so skip it
                                en_d = csp_normalized_slope(sp11, sp22, 0)[0]
                                if en_d == 0:
                                    en_d = csp_normalized_slope(
                                        sp11, sp22, 0.01)[0]
                                if st_d*en_d <= 0:
                                    ray_intersections_count += 1
                                    break
                    else:
                        y1 = csp_at_t(sp1, sp2, t)[1]

                        if y1 == y:
                            # the point is on the path
                            return on_the_path
                        else:
                            if y1 > y and 3*ax*t**2 + 2*bx*t + cx != 0:  # if it's 0 the path only touches the ray
                                ray_intersections_count += 1

    return ray_intersections_count % 2 == 1
