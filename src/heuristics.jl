# (**) Quelques fonctions utilitaires (**)

# Calcule le coût radial d'une route donnée
function route_radius_cost(route::Route, instance::Instance)
    route_diameter = 0
    for id1 in route.order_ids
        for id2 in route.order_ids
            if id1 == id2
                continue
            end
            route_diameter = max(
                route_diameter, instance.euclidean_distances[id1 + 1, id2 + 2]
            )
        end
    end
    vehicle = instance.vehicles[route.family]
    return vehicle.radius_cost * route_diameter / 2
end

# Calcule le coût total d'une route donnée
function route_cost(route::Route, instance::Instance)
    vehicle_idx = route.family
    vehicle = instance.vehicles[vehicle_idx]
    # Coût de location
    cost = vehicle.rental_cost
    # Coût radial
    radius = route_radius_cost(route, instance)
    cost += vehicle.radius_cost * radius
    # Coût de fuel
    n = length(route.order_ids)
    for i in 1:(n-1)
        cost += instance.manhattan_distances[route.order_ids[i] + 1, route.order_ids[i+1] + 1] * vehicle.fuel_cost
    end
    cost += (instance.manhattan_distances[1, route.order_ids[1] + 1] + instance.manhattan_distances[route.order_ids[end] + 1, 1]) * vehicle.fuel_cost
    # Coût total
    return cost
end

# Calcule le poids total d'une route donnée
function route_weight(route::Route, instance::Instance)
    return sum(instance.orders[order_id].weight for order_id in route.order_ids)
end

# Vérifie si une route respecte toutes les contraintes de time windows et capacité
function is_route_feasible(route::Route, instance::Instance)
    vehicle = instance.vehicles[route.family]

    # Vérifier la capacité
    total_weight = route_weight(route, instance)
    if total_weight > vehicle.max_capacity
        return false
    end

    # Vérifier les time windows
    time = 0.0
    current_order_id = instance.depot.id

    for order_id in route.order_ids
        # Temps de trajet
        travel_time = compute_travel_time(current_order_id, order_id, route.family, time, instance)
        time += travel_time

        order = instance.orders[order_id]

        # Attente si on arrive trop tôt
        if time < order.window_start
            time = order.window_start
        end

        # Vérifier qu'on n'arrive pas trop tard
        if time > order.window_end
            return false
        end

        # Temps de service
        time += order.delivery_duration
        current_order_id = order_id
    end

    return true
end

# Teste toutes les orientations possibles de fusion et retourne la meilleure
function best_merge(route1::Route, route2::Route, instance::Instance)
    vehicle_idx = route1.family
    total_weight = route_weight(route1, instance) + route_weight(route2, instance)

    # Vérifier d'abord la capacité
    if total_weight > instance.vehicles[vehicle_idx].max_capacity
        return nothing
    end

    # Tester les 4 orientations possibles
    orientations = [
        (route1.order_ids, route2.order_ids),                    # route1 + route2
        (route1.order_ids, reverse(route2.order_ids)),           # route1 + reverse(route2)
        (reverse(route1.order_ids), route2.order_ids),           # reverse(route1) + route2
        (reverse(route1.order_ids), reverse(route2.order_ids))   # reverse(route1) + reverse(route2)
    ]

    best_cost = Inf
    best_orders = nothing

    for (orders1, orders2) in orientations
        merged_orders = vcat(orders1, orders2)
        candidate_route = Route(vehicle_idx, merged_orders)

        # Vérifier la faisabilité
        if is_route_feasible(candidate_route, instance)
            cost = route_cost(candidate_route, instance)
            if cost < best_cost
                best_cost = cost
                best_orders = merged_orders
            end
        end
    end

    if best_orders === nothing
        return nothing
    end

    return Route(vehicle_idx, best_orders)
end

# Recalcule le meilleur véhicule pour une route donnée (version optimisée)
function optimize_vehicle(route::Route, instance::Instance)
    total_weight = route_weight(route, instance)

    # Pré-filtrer les véhicules par capacité
    candidate_vehicles = [v for v in instance.vehicles
                         if v.max_capacity >= total_weight]

    if isempty(candidate_vehicles)
        return route  # Garder le véhicule actuel si aucun n'est faisable
    end

    best_cost = Inf
    best_family = route.family

    # Tester uniquement les véhicules faisables par capacité
    for vehicle in candidate_vehicles
        candidate_route = Route(vehicle.family, route.order_ids)

        # Vérifier la faisabilité (time windows principalement)
        if is_route_feasible(candidate_route, instance)
            cost = route_cost(candidate_route, instance)
            if cost < best_cost
                best_cost = cost
                best_family = vehicle.family
            end
        end
    end

    return Route(best_family, route.order_ids)
end

# (**) Application de l'algorithme de Clarke-Wright (**)

function clarke_wright_step(solution::Solution, instance::Instance)
    routes = copy(solution.routes)
    nb_routes = length(routes)

    # Matrice des économies avec info sur la route fusionnée
    best_saving = -Inf
    best_i, best_j = -1, -1
    best_merged_route = nothing

    for i in 1:nb_routes
        for j in (i+1):nb_routes
            # Essayer de fusionner les deux routes
            merged_route = best_merge(routes[i], routes[j], instance)

            if merged_route !== nothing
                # Tester aussi avec différents véhicules
                merged_route = optimize_vehicle(merged_route, instance)

                # Calculer l'économie
                cost_before = route_cost(routes[i], instance) + route_cost(routes[j], instance)
                cost_after = route_cost(merged_route, instance)
                saving = cost_before - cost_after

                if saving > best_saving
                    best_saving = saving
                    best_i = i
                    best_j = j
                    best_merged_route = merged_route
                end
            end
        end
    end

    # Si aucune fusion rentable, retourner la solution inchangée
    if best_saving <= 0 || best_merged_route === nothing
        return solution
    end

    # Appliquer la meilleure fusion
    deleteat!(routes, [best_i, best_j])
    push!(routes, best_merged_route)

    return Solution(routes)
end

# (**) Optimisations locales (**)

# Optimisation 2-opt pour une route : inverse un segment pour réduire les croisements
function two_opt_route(route::Route, instance::Instance)
    n = length(route.order_ids)
    if n <= 2
        return route
    end

    improved = true
    best_route = route

    while improved
        improved = false
        current_cost = route_cost(best_route, instance)

        for i in 1:(n-1)
            for j in (i+1):n
                # Créer une nouvelle route en inversant le segment [i, j]
                new_orders = copy(best_route.order_ids)
                new_orders[i:j] = reverse(new_orders[i:j])
                new_route = Route(best_route.family, new_orders)

                # Vérifier la faisabilité et le coût
                if is_route_feasible(new_route, instance)
                    new_cost = route_cost(new_route, instance)
                    if new_cost < current_cost
                        best_route = new_route
                        current_cost = new_cost
                        improved = true
                    end
                end
            end
        end
    end

    return best_route
end

# Applique 2-opt sur toutes les routes d'une solution
function two_opt_solution(solution::Solution, instance::Instance)
    improved_routes = [two_opt_route(route, instance) for route in solution.routes]
    return Solution(improved_routes)
end

# Relocate : déplace un client d'une route vers une autre
function relocate_step(solution::Solution, instance::Instance)
    routes = copy(solution.routes)
    nb_routes = length(routes)

    best_improvement = 0.0
    best_move = nothing

    for i in 1:nb_routes
        route_i = routes[i]
        for pos_i in 1:length(route_i.order_ids)
            customer = route_i.order_ids[pos_i]

            # Essayer d'insérer ce client dans une autre route
            for j in 1:nb_routes
                if i == j
                    continue
                end

                route_j = routes[j]

                # Essayer toutes les positions d'insertion dans route_j
                for pos_j in 0:length(route_j.order_ids)
                    # Nouvelle route i sans le client
                    new_orders_i = [route_i.order_ids[k] for k in 1:length(route_i.order_ids) if k != pos_i]

                    # Nouvelle route j avec le client inséré
                    new_orders_j = copy(route_j.order_ids)
                    insert!(new_orders_j, pos_j + 1, customer)

                    # Créer les nouvelles routes
                    if length(new_orders_i) == 0
                        # Route i devient vide, on la supprime
                        new_route_j = Route(route_j.family, new_orders_j)

                        if is_route_feasible(new_route_j, instance)
                            new_route_j = optimize_vehicle(new_route_j, instance)

                            cost_before = route_cost(route_i, instance) + route_cost(route_j, instance)
                            cost_after = route_cost(new_route_j, instance)
                            improvement = cost_before - cost_after

                            if improvement > best_improvement
                                best_improvement = improvement
                                best_move = (i, j, pos_i, pos_j, new_route_j, nothing)
                            end
                        end
                    else
                        new_route_i = Route(route_i.family, new_orders_i)
                        new_route_j = Route(route_j.family, new_orders_j)

                        if is_route_feasible(new_route_i, instance) && is_route_feasible(new_route_j, instance)
                            new_route_i = optimize_vehicle(new_route_i, instance)
                            new_route_j = optimize_vehicle(new_route_j, instance)

                            cost_before = route_cost(route_i, instance) + route_cost(route_j, instance)
                            cost_after = route_cost(new_route_i, instance) + route_cost(new_route_j, instance)
                            improvement = cost_before - cost_after

                            if improvement > best_improvement
                                best_improvement = improvement
                                best_move = (i, j, pos_i, pos_j, new_route_j, new_route_i)
                            end
                        end
                    end
                end
            end
        end
    end

    # Appliquer le meilleur mouvement
    if best_move !== nothing
        (i, j, pos_i, pos_j, new_route_j, new_route_i) = best_move

        if new_route_i === nothing
            # Route i disparaît
            deleteat!(routes, i)
            # Ajuster l'indice j si nécessaire
            j_adjusted = j > i ? j - 1 : j
            routes[j_adjusted] = new_route_j
        else
            routes[i] = new_route_i
            routes[j] = new_route_j
        end

        return Solution(routes), true
    end

    return solution, false
end

# Exchange : échange deux clients entre deux routes différentes
function exchange_step(solution::Solution, instance::Instance)
    routes = copy(solution.routes)
    nb_routes = length(routes)

    best_improvement = 0.0
    best_move = nothing

    for i in 1:nb_routes
        route_i = routes[i]
        for pos_i in 1:length(route_i.order_ids)
            customer_i = route_i.order_ids[pos_i]

            for j in (i+1):nb_routes
                route_j = routes[j]
                for pos_j in 1:length(route_j.order_ids)
                    customer_j = route_j.order_ids[pos_j]

                    # Échanger les deux clients
                    new_orders_i = copy(route_i.order_ids)
                    new_orders_i[pos_i] = customer_j

                    new_orders_j = copy(route_j.order_ids)
                    new_orders_j[pos_j] = customer_i

                    new_route_i = Route(route_i.family, new_orders_i)
                    new_route_j = Route(route_j.family, new_orders_j)

                    if is_route_feasible(new_route_i, instance) && is_route_feasible(new_route_j, instance)
                        new_route_i = optimize_vehicle(new_route_i, instance)
                        new_route_j = optimize_vehicle(new_route_j, instance)

                        cost_before = route_cost(route_i, instance) + route_cost(route_j, instance)
                        cost_after = route_cost(new_route_i, instance) + route_cost(new_route_j, instance)
                        improvement = cost_before - cost_after

                        if improvement > best_improvement
                            best_improvement = improvement
                            best_move = (i, j, new_route_i, new_route_j)
                        end
                    end
                end
            end
        end
    end

    # Appliquer le meilleur échange
    if best_move !== nothing
        (i, j, new_route_i, new_route_j) = best_move
        routes[i] = new_route_i
        routes[j] = new_route_j
        return Solution(routes), true
    end

    return solution, false
end

# 2-opt* : échange les queues de deux routes
function two_opt_star_step(solution::Solution, instance::Instance)
    routes = copy(solution.routes)
    nb_routes = length(routes)

    best_improvement = 0.0
    best_move = nothing

    for i in 1:nb_routes
        route_i = routes[i]
        for cut_i in 1:(length(route_i.order_ids)-1)
            # Couper route_i en deux parties
            head_i = route_i.order_ids[1:cut_i]
            tail_i = route_i.order_ids[(cut_i+1):end]

            for j in (i+1):nb_routes
                route_j = routes[j]
                for cut_j in 1:(length(route_j.order_ids)-1)
                    # Couper route_j en deux parties
                    head_j = route_j.order_ids[1:cut_j]
                    tail_j = route_j.order_ids[(cut_j+1):end]

                    # Échanger les queues : head_i + tail_j et head_j + tail_i
                    new_orders_i = vcat(head_i, tail_j)
                    new_orders_j = vcat(head_j, tail_i)

                    new_route_i = Route(route_i.family, new_orders_i)
                    new_route_j = Route(route_j.family, new_orders_j)

                    if is_route_feasible(new_route_i, instance) && is_route_feasible(new_route_j, instance)
                        new_route_i = optimize_vehicle(new_route_i, instance)
                        new_route_j = optimize_vehicle(new_route_j, instance)

                        cost_before = route_cost(route_i, instance) + route_cost(route_j, instance)
                        cost_after = route_cost(new_route_i, instance) + route_cost(new_route_j, instance)
                        improvement = cost_before - cost_after

                        if improvement > best_improvement
                            best_improvement = improvement
                            best_move = (i, j, new_route_i, new_route_j)
                        end
                    end
                end
            end
        end
    end

    # Appliquer le meilleur échange
    if best_move !== nothing
        (i, j, new_route_i, new_route_j) = best_move
        routes[i] = new_route_i
        routes[j] = new_route_j
        return Solution(routes), true
    end

    return solution, false
end

# Or-opt : déplace un segment de 1, 2 ou 3 clients
function or_opt_step(solution::Solution, instance::Instance; segment_sizes=[1, 2, 3])
    routes = copy(solution.routes)
    nb_routes = length(routes)

    best_improvement = 0.0
    best_move = nothing

    for i in 1:nb_routes
        route_i = routes[i]

        for segment_size in segment_sizes
            if length(route_i.order_ids) < segment_size
                continue
            end

            # Pour chaque segment de taille segment_size
            for start_pos in 1:(length(route_i.order_ids) - segment_size + 1)
                segment = route_i.order_ids[start_pos:(start_pos + segment_size - 1)]

                # Essayer de le déplacer dans la même route ou une autre
                for j in 1:nb_routes
                    route_j = routes[j]

                    # Positions d'insertion possibles
                    max_insert_pos = (i == j) ? length(route_j.order_ids) - segment_size : length(route_j.order_ids)

                    for insert_pos in 0:max_insert_pos
                        # Ne pas réinsérer au même endroit
                        if i == j && insert_pos >= start_pos - 1 && insert_pos <= start_pos + segment_size - 1
                            continue
                        end

                        # Construire les nouvelles routes
                        if i == j
                            # Mouvement intra-route
                            new_orders = copy(route_i.order_ids)
                            deleteat!(new_orders, start_pos:(start_pos + segment_size - 1))

                            # Ajuster la position d'insertion si nécessaire
                            adjusted_insert_pos = insert_pos >= start_pos ? insert_pos - segment_size : insert_pos

                            for k in 1:segment_size
                                insert!(new_orders, adjusted_insert_pos + k, segment[k])
                            end

                            new_route_i = Route(route_i.family, new_orders)

                            if is_route_feasible(new_route_i, instance)
                                new_route_i = optimize_vehicle(new_route_i, instance)

                                cost_before = route_cost(route_i, instance)
                                cost_after = route_cost(new_route_i, instance)
                                improvement = cost_before - cost_after

                                if improvement > best_improvement
                                    best_improvement = improvement
                                    best_move = (i, j, new_route_i, nothing, true)
                                end
                            end
                        else
                            # Mouvement inter-routes
                            new_orders_i = [route_i.order_ids[k] for k in 1:length(route_i.order_ids)
                                          if k < start_pos || k > start_pos + segment_size - 1]

                            new_orders_j = copy(route_j.order_ids)
                            for k in 1:segment_size
                                insert!(new_orders_j, insert_pos + k, segment[k])
                            end

                            if length(new_orders_i) == 0
                                # Route i disparaît
                                new_route_j = Route(route_j.family, new_orders_j)

                                if is_route_feasible(new_route_j, instance)
                                    new_route_j = optimize_vehicle(new_route_j, instance)

                                    cost_before = route_cost(route_i, instance) + route_cost(route_j, instance)
                                    cost_after = route_cost(new_route_j, instance)
                                    improvement = cost_before - cost_after

                                    if improvement > best_improvement
                                        best_improvement = improvement
                                        best_move = (i, j, nothing, new_route_j, false)
                                    end
                                end
                            else
                                new_route_i = Route(route_i.family, new_orders_i)
                                new_route_j = Route(route_j.family, new_orders_j)

                                if is_route_feasible(new_route_i, instance) && is_route_feasible(new_route_j, instance)
                                    new_route_i = optimize_vehicle(new_route_i, instance)
                                    new_route_j = optimize_vehicle(new_route_j, instance)

                                    cost_before = route_cost(route_i, instance) + route_cost(route_j, instance)
                                    cost_after = route_cost(new_route_i, instance) + route_cost(new_route_j, instance)
                                    improvement = cost_before - cost_after

                                    if improvement > best_improvement
                                        best_improvement = improvement
                                        best_move = (i, j, new_route_i, new_route_j, false)
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # Appliquer le meilleur mouvement
    if best_move !== nothing
        (i, j, new_route_i, new_route_j, intra_route) = best_move

        if intra_route
            routes[i] = new_route_i
        elseif new_route_i === nothing
            # Route i disparaît
            deleteat!(routes, i)
            j_adjusted = j > i ? j - 1 : j
            routes[j_adjusted] = new_route_j
        else
            routes[i] = new_route_i
            routes[j] = new_route_j
        end

        return Solution(routes), true
    end

    return solution, false
end

# Optimisation explicite du radius : réduit les routes trop dispersées
function minimize_radius_step(solution::Solution, instance::Instance)
    routes = copy(solution.routes)

    best_improvement = 0.0
    best_move = nothing

    # Identifier les routes avec un fort radius cost
    for i in 1:length(routes)
        route_i = routes[i]
        vehicle_i = instance.vehicles[route_i.family]

        # Calculer le radius actuel
        max_dist = 0.0
        worst_client = nothing
        worst_client_pos = 0

        for (pos, order_id) in enumerate(route_i.order_ids)
            for other_id in route_i.order_ids
                if order_id != other_id
                    dist = instance.euclidean_distances[order_id + 1, other_id + 1]
                    if dist > max_dist
                        max_dist = dist
                        worst_client = order_id
                        worst_client_pos = pos
                    end
                end
            end
        end

        # Si la route a un fort radius, essayer de déplacer le client le plus éloigné
        radius = max_dist / 2
        radius_penalty = vehicle_i.radius_cost * radius

        # Seuil : optimiser seulement si le radius cost est significatif
        if radius_penalty > 100 && worst_client !== nothing
            # Essayer de déplacer ce client vers une autre route
            for j in 1:length(routes)
                if i == j
                    continue
                end

                route_j = routes[j]

                # Essayer toutes les positions d'insertion
                for insert_pos in 0:length(route_j.order_ids)
                    new_orders_i = [route_i.order_ids[k] for k in 1:length(route_i.order_ids) if k != worst_client_pos]
                    new_orders_j = copy(route_j.order_ids)
                    insert!(new_orders_j, insert_pos + 1, worst_client)

                    if length(new_orders_i) == 0
                        # Route i disparaît
                        new_route_j = Route(route_j.family, new_orders_j)

                        if is_route_feasible(new_route_j, instance)
                            new_route_j = optimize_vehicle(new_route_j, instance)

                            cost_before = route_cost(route_i, instance) + route_cost(route_j, instance)
                            cost_after = route_cost(new_route_j, instance)
                            improvement = cost_before - cost_after

                            if improvement > best_improvement
                                best_improvement = improvement
                                best_move = (i, j, nothing, new_route_j)
                            end
                        end
                    else
                        new_route_i = Route(route_i.family, new_orders_i)
                        new_route_j = Route(route_j.family, new_orders_j)

                        if is_route_feasible(new_route_i, instance) && is_route_feasible(new_route_j, instance)
                            new_route_i = optimize_vehicle(new_route_i, instance)
                            new_route_j = optimize_vehicle(new_route_j, instance)

                            cost_before = route_cost(route_i, instance) + route_cost(route_j, instance)
                            cost_after = route_cost(new_route_i, instance) + route_cost(new_route_j, instance)
                            improvement = cost_before - cost_after

                            if improvement > best_improvement
                                best_improvement = improvement
                                best_move = (i, j, new_route_i, new_route_j)
                            end
                        end
                    end
                end
            end
        end
    end

    # Appliquer le meilleur mouvement
    if best_move !== nothing
        (i, j, new_route_i, new_route_j) = best_move

        if new_route_i === nothing
            deleteat!(routes, i)
            j_adjusted = j > i ? j - 1 : j
            routes[j_adjusted] = new_route_j
        else
            routes[i] = new_route_i
            routes[j] = new_route_j
        end

        return Solution(routes), true
    end

    return solution, false
end

# (**) Heuristiques d'initialisation (**)

# Sweep algorithm : balayage angulaire depuis le dépôt
function sweep_heuristic(instance::Instance)
    routes = Route[]

    # Calculer l'angle de chaque client par rapport au dépôt
    depot = instance.depot
    angles = Float64[]

    for order in instance.orders
        dx = order.longitude - depot.longitude
        dy = order.latitude - depot.latitude
        angle = atan(dy, dx)
        push!(angles, angle)
    end

    # Trier les clients par angle
    sorted_indices = sortperm(angles)

    # Construire des routes en balayant
    current_route_orders = Int[]
    current_vehicle = 1  # Commencer avec le premier véhicule

    for idx in sorted_indices
        order_id = instance.orders[idx].id
        order = instance.orders[order_id]

        # Essayer d'ajouter à la route actuelle
        test_orders = vcat(current_route_orders, [order_id])
        test_route = Route(current_vehicle, test_orders)

        if is_route_feasible(test_route, instance)
            # Ajouter à la route actuelle
            current_route_orders = test_orders
        else
            # Finaliser la route actuelle et en commencer une nouvelle
            if !isempty(current_route_orders)
                final_route = Route(current_vehicle, current_route_orders)
                final_route = optimize_vehicle(final_route, instance)
                push!(routes, final_route)
            end

            # Nouvelle route avec ce client
            current_route_orders = [order_id]

            # Choisir le meilleur véhicule pour ce client
            best_vehicle = 1
            best_cost = Inf
            for vehicle in instance.vehicles
                if vehicle.max_capacity >= order.weight
                    dist = euclidean_distance(depot, order)
                    cost = vehicle.rental_cost + vehicle.fuel_cost * dist
                    if cost < best_cost
                        best_cost = cost
                        best_vehicle = vehicle.family
                    end
                end
            end
            current_vehicle = best_vehicle
        end
    end

    # Finaliser la dernière route
    if !isempty(current_route_orders)
        final_route = Route(current_vehicle, current_route_orders)
        final_route = optimize_vehicle(final_route, instance)
        push!(routes, final_route)
    end

    return Solution(routes)
end

# Regroupement géographique par secteurs
function geographic_clustering(instance::Instance, nb_clusters::Int)
    # Calculer les centroïdes géographiques
    latitudes = [order.latitude for order in instance.orders]
    longitudes = [order.longitude for order in instance.orders]

    lat_min, lat_max = minimum(latitudes), maximum(latitudes)
    lon_min, lon_max = minimum(longitudes), maximum(longitudes)

    # Créer une grille de secteurs
    clusters = [Int[] for _ in 1:nb_clusters]

    # Assigner chaque client au cluster le plus proche (grille simple)
    grid_size = ceil(Int, sqrt(nb_clusters))

    for order in instance.orders
        # Normaliser les coordonnées dans [0, 1]
        lat_norm = (order.latitude - lat_min) / (lat_max - lat_min + 1e-10)
        lon_norm = (order.longitude - lon_min) / (lon_max - lon_min + 1e-10)

        # Trouver la cellule de la grille
        row = min(floor(Int, lat_norm * grid_size) + 1, grid_size)
        col = min(floor(Int, lon_norm * grid_size) + 1, grid_size)

        cluster_id = min((row - 1) * grid_size + col, nb_clusters)
        push!(clusters[cluster_id], order.id)
    end

    # Filtrer les clusters vides
    return filter(!isempty, clusters)
end

# Insertion au plus proche : construit des routes en insérant les clients un par un
function nearest_insertion_heuristic(instance::Instance)
    routes = Route[]
    unvisited = Set(order.id for order in instance.orders)

    # Créer des clusters géographiques
    nb_initial_clusters = max(3, length(instance.orders) ÷ 10)
    clusters = geographic_clustering(instance, nb_initial_clusters)

    # Pour chaque cluster, créer des routes
    for cluster in clusters
        cluster_unvisited = Set(cluster)

        while !isempty(cluster_unvisited)
            # Commencer une nouvelle route avec le client le plus proche du dépôt
            closest_to_depot = nothing
            min_dist = Inf

            for order_id in cluster_unvisited
                order = instance.orders[order_id]
                dist = euclidean_distance(instance.depot, order)
                if dist < min_dist
                    min_dist = dist
                    closest_to_depot = order_id
                end
            end

            # Choisir le meilleur véhicule pour ce client initial
            order = instance.orders[closest_to_depot]
            best_vehicle = nothing
            best_cost = Inf

            for vehicle in instance.vehicles
                if vehicle.max_capacity >= order.weight
                    cost = vehicle.rental_cost + vehicle.fuel_cost * min_dist
                    if cost < best_cost
                        best_cost = cost
                        best_vehicle = vehicle.family
                    end
                end
            end

            current_route = Route(best_vehicle, [closest_to_depot])
            delete!(cluster_unvisited, closest_to_depot)
            delete!(unvisited, closest_to_depot)

            # Insérer les autres clients un par un
            improved = true
            while improved && !isempty(cluster_unvisited)
                improved = false
                best_insertion_cost = Inf
                best_insertion = nothing

                for order_id in cluster_unvisited
                    # Essayer d'insérer à toutes les positions
                    for pos in 0:length(current_route.order_ids)
                        new_orders = copy(current_route.order_ids)
                        insert!(new_orders, pos + 1, order_id)

                        candidate_route = Route(current_route.family, new_orders)

                        if is_route_feasible(candidate_route, instance)
                            cost_increase = route_cost(candidate_route, instance) - route_cost(current_route, instance)

                            if cost_increase < best_insertion_cost
                                best_insertion_cost = cost_increase
                                best_insertion = (order_id, pos, candidate_route)
                                improved = true
                            end
                        end
                    end
                end

                # Appliquer la meilleure insertion
                if best_insertion !== nothing
                    (order_id, pos, candidate_route) = best_insertion
                    current_route = candidate_route
                    delete!(cluster_unvisited, order_id)
                    delete!(unvisited, order_id)
                end
            end

            # Optimiser le véhicule pour cette route
            current_route = optimize_vehicle(current_route, instance)
            push!(routes, current_route)
        end
    end

    # Gérer les clients non visités (hors clusters)
    while !isempty(unvisited)
        order_id = first(unvisited)
        order = instance.orders[order_id]

        best_vehicle = nothing
        best_cost = Inf

        for vehicle in instance.vehicles
            if vehicle.max_capacity >= order.weight
                dist = euclidean_distance(instance.depot, order)
                cost = vehicle.rental_cost + vehicle.fuel_cost * dist
                if cost < best_cost
                    best_cost = cost
                    best_vehicle = vehicle.family
                end
            end
        end

        push!(routes, Route(best_vehicle, [order_id]))
        delete!(unvisited, order_id)
    end

    return Solution(routes)
end

# (**) Variable Neighborhoods Descent (**)

# VND : recherche locale en changeant systématiquement de voisinage
function variable_neighborhood_descent(initial_solution::Solution, instance::Instance)
    solution = initial_solution

    # Liste des voisinages (opérateurs de recherche locale)
    neighborhoods = [
        :relocate,
        :exchange,
        :two_opt_star,
        :or_opt,
        :minimize_radius
    ]

    k = 1  # Indice du voisinage courant

    while k <= length(neighborhoods)
        improved = false

        # Appliquer l'opérateur du voisinage k
        if neighborhoods[k] == :relocate
            new_solution, improved = relocate_step(solution, instance)
        elseif neighborhoods[k] == :exchange
            new_solution, improved = exchange_step(solution, instance)
        elseif neighborhoods[k] == :two_opt_star
            new_solution, improved = two_opt_star_step(solution, instance)
        elseif neighborhoods[k] == :or_opt
            new_solution, improved = or_opt_step(solution, instance)
        elseif neighborhoods[k] == :minimize_radius
            new_solution, improved = minimize_radius_step(solution, instance)
        end

        if improved
            # Amélioration trouvée : accepter et revenir au premier voisinage
            solution = new_solution
            k = 1
        else
            # Pas d'amélioration : passer au voisinage suivant
            k += 1
        end
    end

    return solution
end

# Heuristique avec VND
function vnd_heuristic(instance::Instance)
    # Tester 2 initialisations différentes et garder la meilleure
    solutions = Solution[]

    # Initialisation 1 : nearest insertion avec clustering
    push!(solutions, nearest_insertion_heuristic(instance))

    # Initialisation 2 : sweep algorithm
    push!(solutions, sweep_heuristic(instance))

    # Garder la meilleure initialisation
    best_init_cost = Inf
    solution = solutions[1]
    for sol in solutions
        c = cost(sol, instance)
        if c < best_init_cost
            best_init_cost = c
            solution = sol
        end
    end

    println("    Best initial solution cost: $(round(best_init_cost, digits=2))")

    # Phase 2 : Clarke & Wright
    max_iterations = 100
    iteration = 0
    while iteration < max_iterations
        new_solution = clarke_wright_step(solution, instance)
        if new_solution.routes == solution.routes
            break
        end
        solution = new_solution
        iteration += 1
    end

    # Phase 3 : Optimisations locales (2-opt)
    solution = two_opt_solution(solution, instance)

    # Phase 4 : VND - recherche locale systématique
    # Le VND change automatiquement de voisinage jusqu'à épuisement
    solution = variable_neighborhood_descent(solution, instance)

    # Phase 5 : Dernier passage de 2-opt pour polir
    solution = two_opt_solution(solution, instance)

    return solution
end


