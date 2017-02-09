

# Monter une image de VM sur une partition haswell

pcocc alloc -J -p haswell -c 32 inti-compute 
pcocc alloc -J -p haswell -c 32 profile_dump
# 
pcocc console
pcocc template list
pcocc ssh -l root vm0
pcocc ssh root@vm0
