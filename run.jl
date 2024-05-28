
r = "mp"
instances = ["3", "5", "7"]
timelimit = 3

### ORIGINAL
cd("/Users/vedabradley/Documents/University/RUC/Semester 6/Bachelor Project/Code/Clone/RoRoDualCyclingBreaks/1. Original")
include("1. Original/original.jl"), include("1. Original/original single.jl"), include("1. Original/original grka.jl")

for i in instances
    dualcycling(rule=r, instance=i, timelimit=timelimit, warmstart="single")
end

for i in instances
    dualcycling(rule=r, instance=i, timelimit=timelimit, warmstart="grka", sortmethod="kids")
end

for i in instances
    dualcycling(rule=r, instance=i, timelimit=timelimit, warmstart="grka", sortmethod="top")
end

### SIMPLE BREAK
cd("/Users/vedabradley/Documents/University/RUC/Semester 6/Bachelor Project/Code/Clone/RoRoDualCyclingBreaks/2. Simple Break")
include("2. Simple Break/simple break.jl"), include("2. Simple Break/simple break single.jl"), include("2. Simple Break/simple break grka.jl")

bs = 30
bl = 8

for i in instances
    dualcycling(rule=r, instance=i, timelimit=timelimit, warmstart="grka", sortmethod="kids", breakstart = bs, breaklength = bl)
end

for i in instances
    dualcycling(rule=r, instance=i, timelimit=timelimit, warmstart="grka", sortmethod="top", breakstart = bs, breaklength = bl)
end


### OVERLAPPING BREAK
cd("/Users/vedabradley/Documents/University/RUC/Semester 6/Bachelor Project/Code/Clone/RoRoDualCyclingBreaks/3. Overlapping Break")
include("3. Overlapping Break/overlap break.jl"), include("3. Overlapping Break/overlap break single.jl"), include("3. Overlapping Break/overlap break grka.jl")

bs = 30
bl = 8

for i in instances
    dualcycling(rule=r, instance=i, timelimit=timelimit, warmstart="grka", sortmethod="kids", breakstart = bs, breaklength = bl)
end

for i in instances
    dualcycling(rule=r, instance=i, timelimit=timelimit, warmstart="grka", sortmethod="top", breakstart = bs, breaklength = bl)
end


### TUG ID
cd("/Users/vedabradley/Documents/University/RUC/Semester 6/Bachelor Project/Code/Clone/RoRoDualCyclingBreaks/4. Tug ID")
include("4. Tug ID/tug id.jl"), include("4. Tug ID/tug id single.jl"), include("4. Tug ID/tug id grka.jl")

for i in instances
    dualcycling(rule=r, instance=i, timelimit=timelimit, warmstart="grka", sortmethod="kids")
end

for i in instances
    dualcycling(rule=r, instance=i, timelimit=timelimit, warmstart="grka", sortmethod="top")
end

### LONG BREAK
cd("/Users/vedabradley/Documents/University/RUC/Semester 6/Bachelor Project/Code/Clone/RoRoDualCyclingBreaks/5. Long Break")
include("5. Long Break/long break.jl"), include("5. Long Break/long break single.jl"), include("5. Long Break/long break grka.jl")

bl = 8

for i in instances
    dualcycling(rule=r, instance=i, timelimit=timelimit, warmstart="grka", sortmethod="kids", break_length = bl)
end

for i in instances
    dualcycling(rule=r, instance=i, timelimit=timelimit, warmstart="grka", sortmethod="top", break_length = bl)
end

### SHORT BREAKS
cd("/Users/vedabradley/Documents/University/RUC/Semester 6/Bachelor Project/Code/Clone/RoRoDualCyclingBreaks/6. Short Breaks")
include("6. Short Breaks/short breaks.jl"), include("6. Short Breaks/short breaks single.jl"), include("6. Short Breaks/short breaks grka.jl")

bl = 3
rl = 20
ru = 26

for i in instances
    dualcycling(rule=r, instance=i, timelimit=timelimit, warmstart="grka", sortmethod="kids", breaklength = bl, rl = rl, ru = ru)
end

for i in instances
    dualcycling(rule=r, instance=i, timelimit=timelimit, warmstart="grka", sortmethod="top", breaklength = bl, rl = rl, ru = ru)
end

### COMBINED BREAKS
cd("/Users/vedabradley/Documents/University/RUC/Semester 6/Bachelor Project/Code/Clone/RoRoDualCyclingBreaks/7. Combined Breaks")
include("7. Combined Breaks/combined breaks.jl"), include("7. Combined Breaks/combined breaks single.jl"), include("7. Combined Breaks/combined breaks grka.jl")

sbl = 3
rl = 20
ru = 26

lbl = 8
lbs = 30
lbe = lbs + 14

for i in instances
    dualcycling(rule=r, instance=i, timelimit=timelimit, warmstart="grka", sortmethod="kids", short_breaklength = sbl, rl = rl, ru = ru, long_breaklength = lbl, long_break_period_start = lbs, long_break_period_end = lbe)
end

for i in instances
    dualcycling(rule=r, instance=i, timelimit=timelimit, warmstart="grka", sortmethod="top", short_breaklength = sbl, rl = rl, ru = ru, long_breaklength = lbl, long_break_period_start = lbs, long_break_period_end = lbe)
end


