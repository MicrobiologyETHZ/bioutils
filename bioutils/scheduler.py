#!/usr/bin/env python3
"""
Literature Monitor Scheduler

Automated scheduling system for the literature monitoring tool with email notifications.
"""

import schedule
import time
import subprocess
import smtplib
import os
import json
from datetime import datetime
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email import encoders
import argparse
import sys


class LiteratureScheduler:
    """Automated scheduler for literature monitoring"""
    
    def __init__(self, config_file=None):
        """Initialize scheduler with configuration"""
        # Default configuration
        self.config = {
            "schedule": {
                "frequency": "daily",  # daily, weekly, or custom
                "time": "09:00",       # Time to run (24-hour format)
                "weekday": "monday"    # For weekly runs
            },
            "literature_monitor": {
                "config_file": None,
                "days_back": 1,        # For daily runs, only look at last day
                "min_score": 3,
                "sources": ["pubmed", "biorxiv"],
                "output_dir": "./reports"
            },
            "notifications": {
                "enabled": False,
                "email": {
                    "smtp_server": "smtp.gmail.com",
                    "smtp_port": 587,
                    "sender_email": "",
                    "sender_password": "",  # Use app password for Gmail
                    "recipient_emails": [],
                    "subject_template": "Literature Monitor Report - {date}"
                }
            },
            "retention": {
                "keep_days": 30  # Keep reports for 30 days
            }
        }
        
        # Load custom configuration
        if config_file and os.path.exists(config_file):
            with open(config_file, 'r') as f:
                custom_config = json.load(f)
                self._update_config(self.config, custom_config)
        
        # Ensure output directory exists
        os.makedirs(self.config["literature_monitor"]["output_dir"], exist_ok=True)
    
    def _update_config(self, base_config, update_config):
        """Recursively update configuration"""
        for key, value in update_config.items():
            if key in base_config and isinstance(base_config[key], dict) and isinstance(value, dict):
                self._update_config(base_config[key], value)
            else:
                base_config[key] = value
    
    def run_literature_monitor(self):
        """Execute the literature monitoring script"""
        print(f"[{datetime.now()}] Starting literature monitoring...")
        
        try:
            # Prepare output filename with timestamp
            timestamp = datetime.now().strftime("%Y-%m-%d_%H%M")
            output_file = os.path.join(
                self.config["literature_monitor"]["output_dir"],
                f"literature_report_{timestamp}.csv"
            )
            
            # Build command
            cmd = [
                sys.executable, "-m", "bioutils.literature_monitor",
                "--days", str(self.config["literature_monitor"]["days_back"]),
                "--min-score", str(self.config["literature_monitor"]["min_score"]),
                "--output", output_file,
                "--limit", "50"  # Get more results for scheduled runs
            ]
            
            # Add sources if specified
            if self.config["literature_monitor"]["sources"]:
                cmd.extend(["--sources"] + self.config["literature_monitor"]["sources"])
            
            # Add config file if specified
            if self.config["literature_monitor"]["config_file"]:
                cmd.extend(["--config", self.config["literature_monitor"]["config_file"]])
            
            # Run the literature monitor
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)  # 30 min timeout
            
            if result.returncode == 0:
                print(f"[{datetime.now()}] Literature monitoring completed successfully")
                
                # Send notification if enabled
                if self.config["notifications"]["enabled"]:
                    self._send_notification(output_file, result.stdout)
                
                # Clean up old reports
                self._cleanup_old_reports()
                
                return True
            else:
                print(f"[{datetime.now()}] Literature monitoring failed:")
                print(f"STDOUT: {result.stdout}")
                print(f"STDERR: {result.stderr}")
                return False
                
        except subprocess.TimeoutExpired:
            print(f"[{datetime.now()}] Literature monitoring timed out")
            return False
        except Exception as e:
            print(f"[{datetime.now()}] Error running literature monitoring: {e}")
            return False
    
    def _send_notification(self, output_file, stdout_content):
        """Send email notification with results"""
        if not self.config["notifications"]["email"]["sender_email"]:
            print("Email notifications enabled but no sender email configured")
            return
        
        try:
            # Create message
            msg = MIMEMultipart()
            msg['From'] = self.config["notifications"]["email"]["sender_email"]
            msg['To'] = ", ".join(self.config["notifications"]["email"]["recipient_emails"])
            msg['Subject'] = self.config["notifications"]["email"]["subject_template"].format(
                date=datetime.now().strftime("%Y-%m-%d")
            )
            
            # Email body
            body = f"""
Literature Monitoring Report - {datetime.now().strftime("%Y-%m-%d %H:%M")}

The automated literature monitoring has completed successfully.

Summary from stdout:
{stdout_content[-1000:]}  # Last 1000 characters

Please find the detailed report attached.

Best regards,
Literature Monitor Bot
"""
            
            msg.attach(MIMEText(body, 'plain'))
            
            # Attach CSV files if they exist
            for suffix in ['', '_simple']:
                file_path = output_file.replace('.csv', f'{suffix}.csv')
                if os.path.exists(file_path):
                    with open(file_path, "rb") as attachment:
                        part = MIMEBase('application', 'octet-stream')
                        part.set_payload(attachment.read())
                    
                    encoders.encode_base64(part)
                    part.add_header(
                        'Content-Disposition',
                        f'attachment; filename= {os.path.basename(file_path)}'
                    )
                    msg.attach(part)
            
            # Send email
            server = smtplib.SMTP(
                self.config["notifications"]["email"]["smtp_server"],
                self.config["notifications"]["email"]["smtp_port"]
            )
            server.starttls()
            server.login(
                self.config["notifications"]["email"]["sender_email"],
                self.config["notifications"]["email"]["sender_password"]
            )
            
            text = msg.as_string()
            server.sendmail(
                self.config["notifications"]["email"]["sender_email"],
                self.config["notifications"]["email"]["recipient_emails"],
                text
            )
            server.quit()
            
            print(f"[{datetime.now()}] Email notification sent successfully")
            
        except Exception as e:
            print(f"[{datetime.now()}] Failed to send email notification: {e}")
    
    def _cleanup_old_reports(self):
        """Remove old report files"""
        try:
            output_dir = self.config["literature_monitor"]["output_dir"]
            keep_days = self.config["retention"]["keep_days"]
            cutoff_time = time.time() - (keep_days * 24 * 60 * 60)
            
            for filename in os.listdir(output_dir):
                if filename.startswith("literature_report_") and filename.endswith(".csv"):
                    file_path = os.path.join(output_dir, filename)
                    if os.path.getctime(file_path) < cutoff_time:
                        os.remove(file_path)
                        print(f"[{datetime.now()}] Removed old report: {filename}")
                        
        except Exception as e:
            print(f"[{datetime.now()}] Error cleaning up old reports: {e}")
    
    def setup_schedule(self):
        """Set up the automated schedule"""
        frequency = self.config["schedule"]["frequency"]
        run_time = self.config["schedule"]["time"]
        
        if frequency == "daily":
            schedule.every().day.at(run_time).do(self.run_literature_monitor)
            print(f"Scheduled daily literature monitoring at {run_time}")
            
        elif frequency == "weekly":
            weekday = self.config["schedule"]["weekday"].lower()
            if weekday == "monday":
                schedule.every().monday.at(run_time).do(self.run_literature_monitor)
            elif weekday == "tuesday":
                schedule.every().tuesday.at(run_time).do(self.run_literature_monitor)
            elif weekday == "wednesday":
                schedule.every().wednesday.at(run_time).do(self.run_literature_monitor)
            elif weekday == "thursday":
                schedule.every().thursday.at(run_time).do(self.run_literature_monitor)
            elif weekday == "friday":
                schedule.every().friday.at(run_time).do(self.run_literature_monitor)
            elif weekday == "saturday":
                schedule.every().saturday.at(run_time).do(self.run_literature_monitor)
            elif weekday == "sunday":
                schedule.every().sunday.at(run_time).do(self.run_literature_monitor)
            
            print(f"Scheduled weekly literature monitoring on {weekday}s at {run_time}")
            
            # For weekly runs, look back a full week
            self.config["literature_monitor"]["days_back"] = 7
        
        else:
            print(f"Unknown schedule frequency: {frequency}")
            return False
        
        return True
    
    def run_scheduler(self):
        """Start the scheduler daemon"""
        print(f"[{datetime.now()}] Literature Monitor Scheduler starting...")
        print(f"Next run scheduled for: {schedule.next_run()}")
        
        try:
            while True:
                schedule.run_pending()
                time.sleep(60)  # Check every minute
                
        except KeyboardInterrupt:
            print(f"\n[{datetime.now()}] Scheduler stopped by user")
        except Exception as e:
            print(f"[{datetime.now()}] Scheduler error: {e}")
    
    def run_once(self):
        """Run literature monitoring once (for testing)"""
        return self.run_literature_monitor()
    
    def save_config(self, filename):
        """Save current configuration to file"""
        with open(filename, 'w') as f:
            json.dump(self.config, f, indent=2)
        print(f"Scheduler configuration saved to {filename}")


def create_sample_config():
    """Create a sample configuration file"""
    sample_config = {
        "schedule": {
            "frequency": "daily",
            "time": "09:00",
            "weekday": "monday"
        },
        "literature_monitor": {
            "config_file": "literature_monitor_config.json",
            "days_back": 1,
            "min_score": 3,
            "sources": ["pubmed", "biorxiv"],
            "output_dir": "./reports"
        },
        "notifications": {
            "enabled": True,
            "email": {
                "smtp_server": "smtp.gmail.com",
                "smtp_port": 587,
                "sender_email": "your-email@gmail.com",
                "sender_password": "your-app-password",
                "recipient_emails": ["recipient@email.com"],
                "subject_template": "Literature Monitor Report - {date}"
            }
        },
        "retention": {
            "keep_days": 30
        }
    }
    
    with open("scheduler_config.json", 'w') as f:
        json.dump(sample_config, f, indent=2)
    
    print("Sample scheduler configuration created: scheduler_config.json")
    print("\nIMPORTANT: Edit the configuration file to:")
    print("1. Add your email credentials")
    print("2. Set recipient email addresses")
    print("3. Adjust schedule as needed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Literature Monitor Scheduler')
    parser.add_argument('--config', help='Path to scheduler configuration file')
    parser.add_argument('--run-once', action='store_true', help='Run once instead of scheduling')
    parser.add_argument('--create-config', action='store_true', help='Create sample configuration file')
    
    args = parser.parse_args()
    
    if args.create_config:
        create_sample_config()
        sys.exit(0)
    
    # Create scheduler
    scheduler = LiteratureScheduler(config_file=args.config)
    
    if args.run_once:
        # Run once for testing
        success = scheduler.run_once()
        sys.exit(0 if success else 1)
    else:
        # Set up and run the scheduler
        if scheduler.setup_schedule():
            scheduler.run_scheduler()
        else:
            print("Failed to set up schedule")
            sys.exit(1)