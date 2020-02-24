function sel_turns = findTurns(sol_turns,edges)

not_sel_turns = round(edges(sol_turns <= .49,:))
sel_turns = round(edges(sol_turns >= .99,:))


end

