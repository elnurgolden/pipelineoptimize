import streamlit as st
import numpy as np
from scipy.optimize import linprog
import plotly.graph_objects as go

# Define country names and default values
countries = {
    'Italy': {'price': 12.0, 'capacity': 300.0, 'fuel_cost': 7.0},
    'Bulgaria': {'price': 10.0, 'capacity': 400.0, 'fuel_cost': 3.0},
    'Turkey': {'price': 9.0, 'capacity': 200.0, 'fuel_cost': 2.0},
}

# Function to get input data from user
def get_pipeline_data():
    # Create two columns: left for upstream volume, right for delivery points
    col1, col2 = st.columns([1, 2])  # Adjust the proportions of the columns as needed

    with col1:
        st.write("### Volume Entering the Pipeline")
        # Get the total upstream production volume (entry point)
        upstream_volume = st.number_input("Enter the total upstream production volume (MMBTU):", min_value=0.0, value=800.0, step=0.1)

    with col2:
        st.write("### Delivery Points Information")
        # Initialize lists to store data for each delivery point
        prices = []
        capacities = []
        fuel_costs = []

        # Create input fields for each country
        for country, defaults in countries.items():
            st.write(f"#### {country}")
            price = st.number_input(f"Price of gas at {country} ($/MMBTU):", min_value=0.0, value=defaults['price'], step=0.1, key=f'price_{country}')
            capacity = st.number_input(f"Capacity for {country} (MMBTU):", min_value=0.0, value=defaults['capacity'], step=0.1, key=f'capacity_{country}')
            fuel_cost = st.number_input(f"Fuel (Cost of transportation) for {country} ($/MMBTU):", min_value=0.0, value=defaults['fuel_cost'], step=0.1, key=f'fuel_cost_{country}')
            
            prices.append(price)
            capacities.append(capacity)
            fuel_costs.append(fuel_cost)

    return upstream_volume, prices, capacities, fuel_costs

# Function to calculate non-optimized flow
def calculate_non_optimized_flow(upstream_volume, capacities):
    # Calculate the non-optimized flow considering capacities
    total_capacity = sum(capacities)
    if total_capacity == 0:
        return [0] * len(capacities)
    if upstream_volume <= total_capacity:
        flow = [upstream_volume * (capacity / total_capacity) for capacity in capacities]
    else:
        flow = [min(capacity, upstream_volume * (capacity / total_capacity)) for capacity in capacities]
    return flow

# Function to maximize the value using linear programming
def maximize_value(upstream_volume, prices, capacities, fuel_costs):
    # Objective function: maximize value (price - fuel costs), which is equivalent to minimizing the negative value
    net_prices = [p - fc for p, fc in zip(prices, fuel_costs)]
    c = [-np for np in net_prices]

    # Inequality constraints matrix: add the total upstream production volume constraint (sum of all flows <= upstream_volume)
    A_ub = [[1] * len(capacities)]
    b_ub = [upstream_volume]

    # Bounds for each delivery point (0 <= flow <= capacity)
    bounds = [(0, capacities[i]) for i in range(len(capacities))]

    # Solve the linear program
    result = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method='highs')

    if result.success:
        return result.x, -result.fun, net_prices  # Return the optimized flow, maximum value, and net prices
    else:
        return None, None, None  # In case of failure

# Streamlit app
def main():
    # Set page layout to wide mode
    st.set_page_config(layout="wide")

    st.title("Pipeline Optimization Program")

    # Get pipeline data from the user
    upstream_volume, prices, capacities, fuel_costs = get_pipeline_data()

    # Calculate non-optimized flow
    non_optimized_flow = calculate_non_optimized_flow(upstream_volume, capacities)
    non_opt_total_values = [flow * (price - fc) for flow, price, fc in zip(non_optimized_flow, prices, fuel_costs)]

    # Display non-optimized inputs and results
    st.write("### Non-Optimized Inputs and Results")
    for i, (flow, price, fc) in enumerate(zip(non_optimized_flow, prices, fuel_costs)):
        total_value = flow * (price - fc)
        st.write(f"Delivery Point {i+1} (Non-Optimized):")
        st.write(f"  - Flow: {flow:.2f} MMBTU")
        st.write(f"  - Total Value: ${total_value:.2f}")

    # Maximize the value
    if st.button("Optimize"):
        optimized_flow, max_value, net_prices = maximize_value(upstream_volume, prices, capacities, fuel_costs)

        # Display the results
        if optimized_flow is not None:
            st.success(f"Optimization successful! Maximum Total Value: ${max_value:.2f}")

            # Display optimized flow and total values
            st.write("### Optimized Inputs and Results")
            for i, (flow, net_price) in enumerate(zip(optimized_flow, net_prices)):
                total_value = flow * net_price
                st.write(f"Delivery Point {i+1} (Optimized):")
                st.write(f"  - Flow: {flow:.2f} MMBTU")
                st.write(f"  - Total Value: ${total_value:.2f}")

            # Plot combined chart
            fig = go.Figure()
            # Add non-optimized flow to chart
            fig.add_trace(go.Bar(
                x=[f"Delivery Point {i+1}" for i in range(len(non_optimized_flow))],
                y=non_optimized_flow,
                name='Non-Optimized Flow',
                marker_color='coral'
            ))
            # Add optimized flow to chart
            fig.add_trace(go.Bar(
                x=[f"Delivery Point {i+1}" for i in range(len(optimized_flow))],
                y=optimized_flow,
                name='Optimized Flow',
                marker_color='dodgerblue'
            ))



            fig.update_layout(
                title="Optimized vs. Non-Optimized Flow Distribution",
                xaxis_title="Delivery Points",
                yaxis_title="Flow (MMBTU)",
                barmode='group'
            )
            
            st.plotly_chart(fig)

            # Display total revenues
            st.write("### Total Revenue Comparison")

            # Plot total revenue comparison
            fig2 = go.Figure()
            # Add non-optimized total revenue to chart
            fig2.add_trace(go.Bar(
                x=['Non-Optimized'],
                y=[sum(non_opt_total_values)],
                name='Non-Optimized Revenue',
                marker_color='salmon'
            ))
            # Add optimized total revenue to chart
            fig2.add_trace(go.Bar(
                x=['Optimized'],
                y=[max_value],
                name='Optimized Revenue',
                marker_color='lightgreen'
            ))



            fig2.update_layout(
                title="Total Revenue Comparison",
                xaxis_title="Scenario",
                yaxis_title="Total Revenue ($)",
                barmode='group'
            )
            
            st.plotly_chart(fig2)

        else:
            st.error("Optimization failed. Please check your inputs.")

if __name__ == "__main__":
    main()
