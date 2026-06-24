Get-ChildItem -Recurse -Path . -Exclude out -Include *_genmod.f90, *.mod, *_genmod.smod | Remove-Item -Force
