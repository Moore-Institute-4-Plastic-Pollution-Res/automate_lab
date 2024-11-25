ui <- dashboardPage(
  dashboardHeader(
    title = "Spectral Analysis"
  ),
  
  # Sidebar
  dashboardSidebar(
    sidebarMenu(
      id = "sidebarmenu",
      menuItem(
        "Analyze",
        tabName = "analyze",
        icon = icon("chart-line")
      ),
      menuItem(
        "Instructions",
        tabName = "instructions",
        icon = icon("chalkboard-user")
      )
    )
  ),
  dashboardBody()
)
  