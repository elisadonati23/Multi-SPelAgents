#MWE 

using Agents
using Distributed

@agent Person NoSpaceAgent begin
    type::Symbol
    height::Int64
end

properties = Dict(:time_sim => 0, 
                    :height_to_die => 30, 
                    :height_to_generate => 90,
                    :sum_h => 0.0)



#step functions
function random_h!(Person, model)
    Person.height = rand(1:100)
end

function die!(Person, model)
    if Person.height < model.height_to_die
        remove_agent!(Person, model)
    end
end

function add_people!(Person, model)
    #youngs generate less people
    if Person.type == :young && Person.height > model.height_to_generate
        for i in 1:2
            add_agent!(model, rand((:baby, :young, :adult)), rand(1:100))
        end
    #adults generate more people
    elseif Person.type == :adult && Person.height > model.height_to_generate
        for i in 2:5
            add_agent!(model, rand((:baby, :young, :adult)), rand(1:100))
        end
    end
end

function evolve_environment!(model)
    model.time_sim += 1
    model.sum_h = 0
    for person in allagents(model)
        model.sum_h += getproperty(person, :height)
    end
end

#complex step
function complex_step!(model)

    #first check of agents
    #The collect(values(allagents(model))) is used to create a copy of the agents in the model. 
    #This is necessary because you can't add or remove agents 
    #while iterating over them directly?

    all_agents = collect(values(allagents(model)))
    Threads.@threads for person in all_agents
                        random_h!(person, model)  
                     end
    println(nagents(model))
    
    all_agents = collect(values(allagents(model)))
    young_adults = filter!(person -> person.type == :young || person.type == :adult, all_agents)
    
    #update of agents who survived and can reproduce; only young and adults can generate
    for person in young_adults
        die!(person, model)
        add_people!(person, model)
    end
    
    evolve_environment!(model)
end

model = ABM(Person; properties)

#add agents
for i in 1:100
    add_agent!(model, rand((:baby, :young, :adult)), rand(1:100))
end
#run
step!(model, dummystep, complex_step!, 1000)
allagents(model)