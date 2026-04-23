import numpy as np
from stl import mesh

class Screw:
    def __init__(self, inner_radius, outer_radius, pitch, crest_length, root_length, screw_length, taper=0, tolerance_inset=0):
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.pitch = pitch
        self.crest_length = crest_length
        self.root_length = root_length
        self.screw_length = screw_length
        self.taper = taper
        self.tolerance_inset = tolerance_inset

    def generate_threaded_screw(self):
        """Generate helical threading points for a screw"""
        triangles = []
        num_turns = self.screw_length / self.pitch
        points_per_turn = 12  # Resolution of the helix
        total_points = int(num_turns * points_per_turn)
        
        # Generate points along the helix
        points = []
        for i in range(total_points + 1):
            z = (i / points_per_turn) * self.pitch
            if z > self.screw_length:
                z = self.screw_length
            
            angle = 2 * np.pi * i / points_per_turn
            
            # Create the thread profile with crest and root
            # Crest (outer edge)
            crest_z = z + self.crest_length / 2
            if crest_z <= self.screw_length:
                points.append([
                    self.outer_radius * np.cos(angle),
                    self.outer_radius * np.sin(angle),
                    crest_z
                ])
            
            # Root (inner edge)
            root_z = z - self.root_length / 2
            if root_z >= 0:
                points.append([
                    self.inner_radius * np.cos(angle),
                    self.inner_radius * np.sin(angle),
                    root_z
                ])
        
        return np.array(points)

    def get_mesh(self):
        """Generate mesh from screw points"""
        screw_points = self.generate_threaded_screw()
        
        if len(screw_points) < 3:
            return None
        
        # Create triangles from points (simplified triangulation)
        triangles = []
        num_points = len(screw_points)
        
        for i in range(num_points - 2):
            triangles.append([i, i + 1, i + 2])
        
        # Create mesh object
        screw_mesh = mesh.Mesh(np.zeros(len(triangles), dtype=mesh.Mesh.dtype))
        
        for idx, triangle in enumerate(triangles):
            for i in range(3):
                screw_mesh.vectors[idx][i] = screw_points[triangle[i]]
        
        return screw_mesh


class Cuboid:
    def __init__(self, width, depth, height, screw_inner_radius, screw_outer_radius, taper=0, tolerance_inset=0):
        self.width = width
        self.depth = depth
        self.height = height
        self.screw_inner_radius = screw_inner_radius
        self.screw_outer_radius = screw_outer_radius
        self.taper = taper
        self.tolerance_inset = tolerance_inset

    def generate_cuboid_with_threaded_hole(self):
        """Generate cuboid with threaded hole"""
        # Create cuboid vertices (box centered at origin)
        half_w = self.width / 2
        half_d = self.depth / 2
        
        vertices = np.array([
            # Bottom face
            [-half_w, -half_d, 0],
            [half_w, -half_d, 0],
            [half_w, half_d, 0],
            [-half_w, half_d, 0],
            # Top face
            [-half_w, -half_d, self.height],
            [half_w, -half_d, self.height],
            [half_w, half_d, self.height],
            [-half_w, half_d, self.height],
        ])
        
        # Define faces of the cuboid (2 triangles per face)
        triangles = np.array([
            # Bottom face
            [0, 1, 2], [0, 2, 3],
            # Top face
            [4, 6, 5], [4, 7, 6],
            # Front face
            [0, 5, 1], [0, 4, 5],
            # Back face
            [2, 7, 3], [2, 6, 7],
            # Left face
            [0, 3, 7], [0, 7, 4],
            # Right face
            [1, 5, 6], [1, 6, 2],
        ])
        
        # Add threaded hole representation (cylinder approximation)
        num_segments = 16
        for i in range(num_segments):
            angle1 = 2 * np.pi * i / num_segments
            angle2 = 2 * np.pi * (i + 1) / num_segments
            
            # Top circle points
            p1 = np.array([self.screw_outer_radius * np.cos(angle1), 
                          self.screw_outer_radius * np.sin(angle1), self.height])
            p2 = np.array([self.screw_outer_radius * np.cos(angle2), 
                          self.screw_outer_radius * np.sin(angle2), self.height])
            # Bottom circle points
            p3 = np.array([self.screw_outer_radius * np.cos(angle1), 
                          self.screw_outer_radius * np.sin(angle1), 0])
            p4 = np.array([self.screw_outer_radius * np.cos(angle2), 
                          self.screw_outer_radius * np.sin(angle2), 0])
            
            vertices = np.vstack([vertices, p1, p2, p3, p4])
        
        return vertices, triangles

    def get_mesh(self):
        """Generate mesh from cuboid points"""
        vertices, triangles = self.generate_cuboid_with_threaded_hole()
        
        # Create mesh object
        cuboid_mesh = mesh.Mesh(np.zeros(len(triangles), dtype=mesh.Mesh.dtype))
        
        for idx, triangle in enumerate(triangles):
            for i in range(3):
                if triangle[i] < len(vertices):
                    cuboid_mesh.vectors[idx][i] = vertices[triangle[i]]
        
        return cuboid_mesh


# You can create instances of Screw and Cuboid classes and generate STL files.

if __name__ == "__main__":
    # Example usage
    screw = Screw(inner_radius=5, outer_radius=10, pitch=2, crest_length=3, root_length=1, screw_length=20)
    cuboid = Cuboid(width=15, depth=15, height=10, screw_inner_radius=screw.inner_radius, screw_outer_radius=screw.outer_radius)
    screw_mesh = screw.get_mesh()
    cuboid_mesh = cuboid.get_mesh()
    
    # Save meshes as STL files
    if screw_mesh:
        screw_mesh.save('screw.stl')
        print("Screw mesh saved as screw.stl")
    
    if cuboid_mesh:
        cuboid_mesh.save('cuboid.stl')
        print("Cuboid mesh saved as cuboid.stl")