library(emayili)

run_pipeline_emayili <- function() {
  on.exit({
    cat("📧 Sending via emayili...\n")
    
    # Evaluate variables FIRST
    finish_time <- Sys.time()
    hostname <- Sys.info()[["nodename"]]
    
    email <- envelope() %>%
      from("guruguhan23@gmail.com") %>%
      to("guruguhan23@gmail.com") %>%
      emayili::subject("✅ R Job Done") %>%
      text(sprintf(
        "Pipeline finished at %s\nHost: %s\nReady for review!",
        format(finish_time, "%Y-%m-%d %H:%M:%S"),
        hostname
      ))
    
    smtp <- server(
      host     = "smtp.gmail.com",
      port     = 587,
      username = "guruguhan23@gmail.com",
      password = Sys.getenv("GMAIL_APP_PASSWORD")
    )
    
    smtp(email, verbose = TRUE)
    cat("✅ Email sent with real values!\n")
  }, add = TRUE)
  
  # Your simulation/real job
  cat("Starting job...\n")
  for (i in 1:10) {
    Sys.sleep(1)
    cat(sprintf("Step %d/10\n", i))
  }
  cat("Job completed.\n")
}

run_pipeline_emayili()


# === for html emails ===
email <- envelope() %>%
  from("guruguhan23@gmail.com") %>%
  to("guruguhan23@gmail.com") %>%
  subject("🚀 R Analysis Complete") %>%
  html(sprintf(
    "<h2>✅ Job Finished!</h2>
    <p><strong>Time:</strong> %s</p>
    <p><strong>Host:</strong> %s</p>
    <ul>
      <li>Steps: 10/10</li>
      <li>Status: Success</li>
    </ul>",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    Sys.info()[["nodename"]]
  ))