import numpy as np
import gmsh


def mesh_func(X,
         mesh_size_at_airfoil,
         mesh_size_at_farfield,
         farfield_factor,
         model_name = 'airfoil',
         mesh_file_location = None,
         visualize = False
    ):
    
    # Initialize gmsh and add a model
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 1)
    
    gmsh.model.add(model_name)


    # Initialize list of points on shape and add points
    points_on_shape = []

    total_points_on_shape = X.shape[0] - 1
    point_tag = 1
    for i in range(total_points_on_shape):
        x, y = X[i, 0].item(), X[i, 1].item()

        points_on_shape.append(gmsh.model.geo.addPoint(x, y, 0, mesh_size_at_airfoil, point_tag))

        point_tag += 1
    

    # Initialize list of lines on shape and add lines
    lines_on_shape = []

    total_lines_on_shape = total_points_on_shape
    line_tag = 1
    for i in range(total_lines_on_shape):
        if i == total_lines_on_shape - 1:
            pt1, pt2 = points_on_shape[i], points_on_shape[0]
        else:
            pt1, pt2 = points_on_shape[i], points_on_shape[i + 1]
        
        lines_on_shape.append(gmsh.model.geo.addLine(pt1, pt2, line_tag))

        line_tag += 1
    

    # Create a curve loop for the shape
    curve_loop_tag = 1
    curve_loop_shape = gmsh.model.geo.addCurveLoop(lines_on_shape, curve_loop_tag)
    curve_loop_tag += 1


    # Create far field boundary
    xmin, xmax = -1 * farfield_factor, 1 * farfield_factor
    ymin, ymax = -1 * farfield_factor, 1 * farfield_factor

    # Add boundary points
    points_on_boundary = []
    points_on_boundary.append(gmsh.model.geo.addPoint(xmin, ymin, 0, mesh_size_at_farfield, point_tag))
    point_tag += 1
    points_on_boundary.append(gmsh.model.geo.addPoint(xmax, ymin, 0, mesh_size_at_farfield, point_tag))
    point_tag += 1
    points_on_boundary.append(gmsh.model.geo.addPoint(xmax, ymax, 0, mesh_size_at_farfield, point_tag))
    point_tag += 1
    points_on_boundary.append(gmsh.model.geo.addPoint(xmin, ymax, 0, mesh_size_at_farfield, point_tag))
    point_tag += 1

    # Add boundary lines
    lines_on_boundary = []
    lines_on_boundary.append(gmsh.model.geo.addLine(points_on_boundary[0], points_on_boundary[1], line_tag))
    line_tag += 1
    lines_on_boundary.append(gmsh.model.geo.addLine(points_on_boundary[1], points_on_boundary[2], line_tag))
    line_tag += 1
    lines_on_boundary.append(gmsh.model.geo.addLine(points_on_boundary[2], points_on_boundary[3], line_tag))
    line_tag += 1
    lines_on_boundary.append(gmsh.model.geo.addLine(points_on_boundary[3], points_on_boundary[0], line_tag))
    line_tag += 1


    # Create a curve loop for the boundary
    curve_loop_boundary = gmsh.model.geo.addCurveLoop(lines_on_boundary, curve_loop_tag)
    curve_loop_tag += 1


    # Create a plane surface for the domain
    plane_surface_tag = 1
    plane_surface = gmsh.model.geo.addPlaneSurface([curve_loop_shape, curve_loop_boundary], plane_surface_tag)
    plane_surface_tag += 1


    gmsh.model.geo.synchronize()

    airfoil_curves = lines_on_shape
    farfield_curves = lines_on_boundary

    gmsh.model.geo.synchronize()

    # Add physical groups
    gmsh.model.addPhysicalGroup(1, airfoil_curves, name = 'airfoil')
    gmsh.model.addPhysicalGroup(1, farfield_curves, name = 'farfield')
    gmsh.model.addPhysicalGroup(2, [plane_surface], name = 'plane_surface')

    # Save tags of entities in physical groups in a dictionary
    physical_groups = {
        'airfoil_curves': airfoil_curves,
        'farfield_curves': farfield_curves,
        'plane_surface': plane_surface
    }


    # Generate the mesh
    gmsh.model.mesh.generate()

    # Get total nodes and elements
    nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()
    num_nodes = len(nodeTags)
    _, elementTags, _ = gmsh.model.mesh.getElements()
    num_elements = sum(len(tags) for tags in elementTags)

    mesh_details = {'num_nodes': num_nodes, 'num_elements': num_elements}


    if mesh_file_location is not None:
        gmsh.write(mesh_file_location)

    if visualize:
        gmsh.fltk.run()

    gmsh.finalize()

    return mesh_details, physical_groups