"""
@ Simulations for leapfrogging vortex rings using the Biot-Savart law and Lamb-Ossen vortex core model.
@ author ZimoJupiter
@ w.zimo@outlook.com
@ date 06 Aug 2025
@ license MIT License
"""
import numpy as np
from numpy import cos, sin, tan, pi, sqrt, exp, log
from numpy import arctan, arccos, arctan2, deg2rad, rad2deg
import pandas as pd
from numba import jit
from numba import njit
import tqdm
import matplotlib.pyplot as plt
import vtk
import scipy
import imageio
import glob

plt.rcParams['font.serif'] = 'Times New Roman'
# plt.rcParams['font.weight'] = 'normal'
plt.rcParams["font.size"] = 8
plt.rcParams["figure.figsize"] = (3.25, 3.25*3/4)
plt.rcParams["figure.dpi"] = 600
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": "Times New Roman",
    "text.latex.preamble": r"\usepackage{mathptmx}", })

def main():
    t = np.arange(0, 120.01, 1)
    dt = t[1] - t[0]

    theta = np.arange(0, 2*pi + 2*pi/144, 2*pi/144)
    x = sin(theta)
    y = cos(theta)
    z = 0*np.ones_like(theta)

    Start = np.zeros((2, t.shape[0], theta.shape[0], 3))
    End = np.zeros((2, t.shape[0], theta.shape[0], 3))
    Start[0, :, :, 0] = 0.99*sin(theta)[np.newaxis, :]
    Start[0, :, :, 1] = 0.99*cos(theta)[np.newaxis, :]
    Start[1, :, :, 0] = sin(theta)[np.newaxis, :]
    Start[1, :, :, 1] = cos(theta)[np.newaxis, :]
    Start[0, :, :, 2] = 0
    Start[1, :, :, 2] = 0.2
    End[0] = np.roll(Start[0], shift=-1, axis=1)
    End[1] = np.roll(Start[1], shift=-1, axis=1)

    Gamma = 0.1*np.ones((2, t.shape[0], theta.shape[0]))
    rc = 0.001*np.ones((2, t.shape[0], theta.shape[0]))
    Direction = np.array([1, 1])

    Position = Start.copy()

    # iteration (explicit Euler)
    Vind = np.zeros((2, t.shape[0], theta.shape[0], 3))

    for t_i in range(1, t.shape[0]-1):
        for i in range(2):
            for theta_i in range(theta.shape[0]):
                Vind[i, t_i, theta_i, :] = BiotSavart(Start[:,t_i, :, :], 
                                                        End[:,t_i, :, :],
                                                        Gamma[:,t_i, :], 
                                                        Position[i, t_i, theta_i, :], 
                                                        rc[:,t_i, :], 
                                                        Direction,
                                                        2)

        Position[:, t_i+1, :, :] = Position[:, t_i, :, :] + Vind[:, t_i, :, :]*dt

        Visualization(2, Position[:, t_i, :, :], 100, 
                      f'Plots/3D_Trajectories_Visualization_{round(t[t_i], 2)}.png',
                      time_text=f"time={round(t[t_i], 2)}s")

    file_list = sorted(glob.glob('Plots/3D_Trajectories_Visualization_*.png'))
    file_list = sorted(file_list, key=lambda x: int(x.split('_')[-1].split('.')[0]))
    gif_images = []
    for file in file_list:
        image = imageio.imread(file)
        gif_images.append(image)
    imageio.mimsave(f'Plots/Animation.gif', gif_images, loop=0, duration=10)
    print('Animate is READY!')

    breakpoint

def BiotSavart(Start, End, Gamma, Location, rc, Direction, n=2): 
    # Left hand rule for Biot-Savart law
    # Start, End, Location, Direction: (N, 3) arrays
    # n = 2 # Vortex core: 1,3_Vatistas; 2_Lamb-Ossen; \infinity_Rankine
    R1 = Location - End
    R2 = Location - Start
    R3 = End - Start
    Direction = Direction[(...,) + (np.newaxis,)*(R3.ndim-1)]
    R3 *= Direction
    Vind = np.zeros((Start.shape))

    # Tensors for Lamb-Ossen
    Gamma = np.expand_dims(Gamma, axis=-1)
    rc = np.expand_dims(rc, axis=-1)
    distance = np.linalg.norm(np.cross(R2, R1), axis=-1)/np.linalg.norm(R3, axis=-1)
    cosR1R3 = np.sum(R1*R3, axis=-1) / (np.linalg.norm(R1, axis=-1)*np.linalg.norm(R3, axis=-1))
    cosR2R3 = np.sum(R2*R3, axis=-1) / (np.linalg.norm(R2, axis=-1)*np.linalg.norm(R3, axis=-1))
    crossR2R1 = np.cross(R2,R1)
    normcrossR2R1 = np.linalg.norm(crossR2R1, axis=-1)
    distance = np.expand_dims(distance, axis=-1) # Dimension of the tensors
    cosR1R3 = np.expand_dims(cosR1R3, axis=-1)
    cosR2R3 = np.expand_dims(cosR2R3, axis=-1)
    normcrossR2R1 = np.expand_dims(normcrossR2R1, axis=-1)
    if n == 0:
        h = distance + 1e-6
    else:
        h = (rc**(2*n) + distance**(2*n))**(1/n)/(distance) # Lamb-Ossen, n = 2
    Vind = Gamma/(4*pi*h)*(cosR1R3-cosR2R3)*crossR2R1/normcrossR2R1
    
    Vind[np.isnan(Vind)] = 0 # Clearing points falling on vectors
    Vind = np.sum(Vind, axis=tuple(range(Vind.ndim-1)))
    
    return Vind

def Visualization(N, Rvec_filament, factor, file_name, time_text=None):
    # Visualize the vortex ring trajectories using VTK
    # N: number of vortex rings
    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindow.OffScreenRenderingOn()
    renderer.SetBackground(1.0, 1.0, 1.0)
    colors = np.array([[1, 0, 0], 
                       [0, 0, 1]])

    for bn in range(N):
        x = Rvec_filament[bn, :, 0]
        y = Rvec_filament[bn, :, 1]
        z = Rvec_filament[bn, :, 2]
        points = vtk.vtkPoints()
        for i in range(len(x)):
            points.InsertNextPoint(x[i], y[i], z[i])

        line = vtk.vtkPolyLine()
        line.GetPointIds().SetNumberOfIds(len(x))
        for i in range(len(x)):
            line.GetPointIds().SetId(i, i)
        lines = vtk.vtkCellArray()
        lines.InsertNextCell(line)
        polyData = vtk.vtkPolyData()
        polyData.SetPoints(points)
        polyData.SetLines(lines)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(polyData)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(colors[bn,0], colors[bn,1], colors[bn,2])
        actor.GetProperty().SetLineWidth(0.1*factor)
        renderer.AddActor(actor)

        marker = vtk.vtkSphereSource()
        marker.SetCenter(x[0], y[0], z[0])
        marker.SetRadius(0.005)
        markerMapper = vtk.vtkPolyDataMapper()
        markerMapper.SetInputConnection(marker.GetOutputPort())
        markerActor = vtk.vtkActor()
        markerActor.SetMapper(markerMapper)
        markerActor.GetProperty().SetColor(1, 0, 0)
        renderer.AddActor(markerActor)
    
    if time_text is not None:
        text_actor = vtk.vtkTextActor()
        text_actor.SetInput(time_text)
        text_actor.GetTextProperty().SetColor(0, 0, 0)
        text_actor.GetTextProperty().SetFontSize(20)
        text_actor.GetTextProperty().SetJustificationToRight()
        text_actor.GetTextProperty().SetVerticalJustificationToTop()
        text_actor.SetPosition(550, 430)
        renderer.AddActor2D(text_actor)

    camera = renderer.GetActiveCamera()
    camera.SetViewUp(0, 0, 1)
    camera.SetPosition(0, -10, 2)
    camera.SetFocalPoint(0, 0.5*(max(y)+min(y)), 0.5*(max(z)+min(z)))

    renderWindow.SetSize(600, 450)
    renderWindow.Render()
    windowToImageFilter = vtk.vtkWindowToImageFilter()
    windowToImageFilter.SetInput(renderWindow)
    windowToImageFilter.Update()
    writer = vtk.vtkPNGWriter()
    writer.SetFileName(file_name)
    writer.SetInputConnection(windowToImageFilter.GetOutputPort())
    writer.Write()

if __name__ == "__main__":
    main()