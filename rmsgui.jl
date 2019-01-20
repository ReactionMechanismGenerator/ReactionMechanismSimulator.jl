using Gtk

win1 = GtkWindow("Simulation Properties")
g = GtkGrid()
rlabel = GtkLabel("Reactor Type")
tlabel = GtkLabel("Final Time (sec)")
Tlabel = GtkLabel("T (K)")
Plabel = GtkLabel("P (Pa)")
Vlabel = GtkLabel("V (m^3)")
rtollabel = GtkLabel("rtol")
atollabel = GtkLabel("atol")
nlabels = GtkLabel("Name")
nlabel = GtkLabel("Initial Moles")
Reacttype = GtkComboBoxText()
rchoices = ["ConstantTPDomain","ConstantTVDomain","ConstantVDomain"]
for choice in rchoices
    push!(Reacttype,choice)
end
set_gtk_property(Reacttype,:active,1)
tf = GtkEntry()
T0 = GtkEntry()
P0 = GtkEntry()
V0 = GtkEntry()
rtolent = GtkEntry()
atolent = GtkEntry()
nsNames = Array([GtkEntry() for i = 1:15])
ns = Array([GtkEntry() for i = 1:15])
sens = GtkCheckButton("Sensitivity")
start = GtkButton("Simulate")
g[1,1] = rlabel
g[1,2] = tlabel
g[1,3] = Tlabel
g[1,4] = Plabel
g[1,5] = Vlabel
g[2,1] = Reacttype
g[2,2] = tf
g[2,3] = T0
g[2,4] = P0
g[2,5] = V0
g[3,1] = sens
g[3,2] = start
g[1,6] = rtollabel
g[2,6] = rtolent
g[1,7] = atollabel
g[2,7] = atolent
g[1,8] = nlabels
g[2,8] = nlabel
ns = Array{GtkEntry,1}([])
nsNames = Array{GtkEntry,1}([])
for k = 1:15
    name = GtkEntry()
    n = GtkEntry()
    g[1,8+k] = name
    g[2,8+k] = n
    push!(ns,n)
    push!(nsNames,name)
end

function on_button_clicked(w)
  T = parse(Float64,get_gtk_property(T0,:text,String))
  println(T)
end
signal_connect(on_button_clicked, start, "clicked")

set_gtk_property!(g, :column_homogeneous, true)
set_gtk_property!(g, :column_spacing, 15)
set_gtk_property!(g, :row_spacing, 8)

push!(win1,g)
showall(win1)
