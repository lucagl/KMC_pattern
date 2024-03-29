import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


# L=10
# total_frames =10

conf_file = open ("configuration.txt")
L, total_frames, final_T, final_conc = conf_file.readline().split()
L = int(L)
total_frames = int(total_frames)

print ("Simulation box size =", L)
print ("Number of frames =", total_frames)




curr_pos =0

file = open("island_conv0_sigma2.txt")
dummy,dummy, T,dummy,c  = file.readline().split()
file.close()


print ("initial T = ", T)
print ("initial c=", c)

island = np.loadtxt("island_conv0_sigma2.txt",delimiter='\t')
adatom = np.loadtxt("adatom_conv0_sigma6.txt",delimiter='\t')
z_island = island

Z_island=z_island.reshape((L,L))


#z_attachable = adatom[:,1]


z_adatom = adatom

Z_adatom = z_adatom.reshape(L,L)






def key_event(e):

    global curr_pos
    if (e.key == "right") or (e.key == "up"):
            curr_pos = curr_pos+1
    #print 1
    elif (e.key == "left") or (e.key == "down"):
            curr_pos = curr_pos - 1
    elif (e.key == "d"):
            curr_pos = curr_pos+10
    elif (e.key == "a"):
        curr_pos = curr_pos-10
    else:
        return
    curr_pos = curr_pos % (total_frames+1)
    #print(curr_pos)

    file = open("island_conv" + str(curr_pos) + "_sigma2.txt")
    dummy,dummy, T, dummy,c = file.readline().split()

    file_island = "island_conv" + str(curr_pos) + "_sigma2.txt"
    file_adatom = "adatom_conv" + str(curr_pos) + "_sigma6.txt"

    island = np.loadtxt(file_island,delimiter='\t')
    adatom = np.loadtxt(file_adatom,delimiter='\t')

    z_island = island


    Z_island=z_island.reshape((L,L))


    z_adatom = adatom
    # max= np.amax(z_adatom)
    # print (max)

    Z_adatom = z_adatom.reshape(L,L)

    plt.ion()


    fig.suptitle('frame ' + str(curr_pos) + ' of ' + str(total_frames) + "\n" + "T =" + T + "average conc "+ c)

    island_plot = ax[0].imshow(Z_island, interpolation='none')
    ax[0].set_title('Island')
    ax[0].axis('off')


    adatom_plot = ax[1].imshow(Z_adatom, cmap = 'afmhot' , interpolation = 'none')

    #cbar = plt.colorbar(adatom_plot, cax =cax, orientation = 'horizontal')
    # cbar.set_clim(0, max)

    ax[1].set_title('Adatom density')
    ax[1].axis('off')








fig,ax = plt.subplots(1,2)
fig.suptitle('frame ' + str(curr_pos) + ' of ' + str(total_frames)+ "\n" + "T =" + T + "average conc " + c )

island_plot = ax[0].imshow(Z_island, interpolation='none')
ax[0].set_title('Island')
ax[0].axis('off')
adatom_plot = ax[1].imshow(Z_adatom, cmap ='afmhot' , interpolation = 'none')
ax[1].set_title('Adatom density')
ax[1].axis('off')



#divider = make_axes_locatable(ax[1])
#cax=divider.append_axes("bottom", size='5%', pad = 0.1)

#cbar = plt.colorbar(adatom_plot, cax =cax, orientation = 'horizontal')
#cbar = plt.colorbar(adatom_plot, ax =ax[1], orientation = 'horizontal', fraction = 0.05 )

fig.canvas.mpl_connect('key_press_event', key_event)
plt.show()



# Estetics ---------------------------------------


#fig.colorbar(adatom_plot,ax=ax[1], fraction = 0.05,orientation = 'horizontal')


